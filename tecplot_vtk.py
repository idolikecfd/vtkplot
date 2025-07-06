import vtk
import numpy as np
import argparse
import os.path
import sys

def process_tecplot_file(input_file):
    """
    Process a Tecplot data file and visualize it using VTK.
    
    Parameters:
    -----------
    input_file : str
        Path to the Tecplot .dat file
    """
    # Check if the input file exists
    if not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        sys.exit(1)
        
    print(f"Processing file: {input_file}")
    
    # 1. Load the Tecplot data file
    reader = vtk.vtkTecplotReader()
    reader.SetFileName(input_file)
    reader.Update()

    # Initialize actors list
    actors = []

    # Get the multi-block output
    multiblock_data = reader.GetOutput()
    print(f"Data type: {multiblock_data.GetClassName()}")
    print(f"Number of blocks: {multiblock_data.GetNumberOfBlocks()}")

    # 2. Process each block in the multi-block dataset
    for i in range(multiblock_data.GetNumberOfBlocks()):
        block = multiblock_data.GetBlock(i)
        if block is None:
            print(f"Block {i} is None, skipping")
            continue
        
        print(f"Block {i} type: {block.GetClassName()}")
        
        # Apply contour settings from the Tecplot file
        # In Tecplot we have variable 5 (u) set for contour
        contour = vtk.vtkContourFilter()
        contour.SetInputData(block)
        
        # Set up contour values based on Tecplot levels
        # These match the raw data values in the Tecplot layout
        contour_values = [
            0.335123736549, 0.344903259539, 0.35468278253, 0.364462305521,
            0.374241828511, 0.384021351502, 0.393800874493, 0.403580397483,
            0.413359920474, 0.423139443464, 0.432918966455, 0.442698489446,
            0.452478012436, 0.462257535427, 0.472037058418, 0.481816581408,
            0.491596104399, 0.50137562739, 0.51115515038
        ]
        
        # In Tecplot, VAR = 5 corresponds to 'u' velocity
        # In VTK, array indices start at 0, so we need index 4
        contour.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "u")
        
        for value in contour_values:
            contour.SetValue(contour_values.index(value), value)
        
        contour.Update()
        
        # Create a mapper for the contour output
        contour_mapper = vtk.vtkPolyDataMapper()
        contour_mapper.SetInputConnection(contour.GetOutputPort())
        contour_mapper.SetScalarRange(0.325344213558, 0.520934673371)  # From CMIN/CMAX in Tecplot layout
        
        # Create a lookup table to match Tecplot's coloring as seen in reference image
        lut = vtk.vtkLookupTable()
        # Based on the reference image, we need a color map that goes from 
        # blue -> cyan -> green -> yellow -> orange -> red
        lut.SetNumberOfColors(256)
        # Reverse the hue range to match the Tecplot colors (red for high values, blue for low)
        lut.SetHueRange(0.667, 0.0)  
        lut.SetSaturationRange(1.0, 1.0)
        lut.SetValueRange(1.0, 1.0)
        lut.Build()
        contour_mapper.SetLookupTable(lut)
        
        # Create a contour actor
        contour_actor = vtk.vtkActor()
        contour_actor.SetMapper(contour_mapper)
        actors.append(contour_actor)
        
        # Create a mapper for the original dataset (for mesh edges)
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(block)
        
        # Create an actor for edges only
        edge_actor = vtk.vtkActor()
        edge_actor.SetMapper(mapper)
        
        # Set edge properties according to the Tecplot file
        # EDGELAYER SHOW = YES, COLOR = BLACK, LINETHICKNESS = 0.1
        edge_actor.GetProperty().SetRepresentationToWireframe()
        edge_actor.GetProperty().SetColor(0, 0, 0)  # Black
        edge_actor.GetProperty().SetLineWidth(0.1)
        actors.append(edge_actor)
        
        # Configure surface properties from the Tecplot file
        # Set lighting effects to Gouraud
        contour_actor.GetProperty().SetInterpolationToGouraud()
        
        # Create surface mapper for flood contour
        surface_mapper = vtk.vtkDataSetMapper()
        surface_mapper.SetInputData(block)
        surface_mapper.SetScalarRange(0.325344213558, 0.520934673371)  # Same as contour range
        surface_mapper.SetScalarModeToUsePointFieldData()
        surface_mapper.SelectColorArray("u")  # Variable 5 - 'u' velocity
        surface_mapper.SetLookupTable(lut)
        
        # Create surface actor
        surface_actor = vtk.vtkActor()
        surface_actor.SetMapper(surface_mapper)
        surface_actor.GetProperty().SetInterpolationToGouraud()  # LIGHTINGEFFECT = GOURAUD
        actors.append(surface_actor)

    # 3. Set up the renderer and add all actors
    renderer = vtk.vtkRenderer()
    for actor in actors:
        renderer.AddActor(actor)
    renderer.SetBackground(0.6, 0.3, 0)  # Brown/orange background to match reference image

    # 7. Configure lighting similar to Tecplot settings
    #light = vtk.vtkLight()
    #light.SetPosition(-0.2, -0.2, 0.959)  # From LIGHTSOURCE XYZDIRECTION
    #light.SetIntensity(0.75)  # From INTENSITY = 75
    #renderer.AddLight(light)


    # 4. Set up the render window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(800, 800)  # Matching the layout's WIDTH = 8, HEIGHT = 8 but in pixels

    # 5. Set up camera using values from the Tecplot file
    camera = renderer.GetActiveCamera()
    # Using the values from $!THREEDVIEW and $!GLOBALTHREED sections
    # Set camera position based on VIEWERPOSITION and angles
    # Convert Tecplot's PSI and THETA angles to camera position
    psi_angle = 65.8537   # From PSIANGLE in layout
    theta_angle = 140.921 # From THETAANGLE in layout

    # Adjust these values to match the view in the reference image
    # The position values need adjustment to match the reference view
    camera.SetPosition(-872.35, 1075.64, 620.81)
    camera.SetFocalPoint(0, 0, 0)  # Center on the wing
    camera.SetViewUp(0, 0, 1)      # Z-up orientation as shown in reference
    camera.SetViewAngle(30)        # Adjusted for better view

    # Fine-tune to match the specific wing orientation in the reference image
    camera.Elevation(-20)
    camera.Azimuth(10)

    # 6. Set up the interactor
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Add key bindings to toggle display modes
    style = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(style)

    # Create a color bar / scalar bar for the contour values
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(lut)
    scalar_bar.SetTitle("u velocity")
    scalar_bar.SetNumberOfLabels(10)
    # Position the scalar bar similar to Tecplot's legend position
    scalar_bar.SetPosition(0.8, 0.1)
    scalar_bar.SetPosition2(0.15, 0.8)
    renderer.AddActor2D(scalar_bar)

    # 7. Initialize the interactor and start rendering
    renderWindowInteractor.Initialize()
    renderer.ResetCamera()  # Make sure the camera shows all actors

    # Add axis display to match the reference image (shows X, Y, Z axes in top right)
    axes = vtk.vtkAxesActor()
    axes.SetTotalLength(30, 30, 30)
    axes.SetShaftTypeToLine()
    axes.SetNormalizedShaftLength(1, 1, 1)
    axes.SetNormalizedTipLength(0.1, 0.1, 0.1)
    # Make the axis labels visible
    axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(12)
    axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(12)
    axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(12)

    # Create the orientation widget to display the axes
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOutlineColor(0.9300, 0.5700, 0.1300)
    axes_widget.SetOrientationMarker(axes)
    axes_widget.SetInteractor(renderWindowInteractor)
    axes_widget.SetViewport(0.8, 0.8, 1.0, 1.0)  # Top right corner to match reference
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    renderWindow.Render()
    renderWindowInteractor.Start()

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Convert and visualize Tecplot data files using VTK')
    
    # Add arguments
    parser.add_argument('input_file', type=str, help='Path to the input Tecplot .dat file')
    
    # Optional arguments can be added here
    parser.add_argument('--window-size', type=int, nargs=2, metavar=('WIDTH', 'HEIGHT'),
                        default=[800, 800], help='Window size in pixels (default: 800 800)')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Process the Tecplot file
    process_tecplot_file(args.input_file)

if __name__ == '__main__':
    main()
