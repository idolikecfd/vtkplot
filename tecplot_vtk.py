import vtk
import numpy as np
import argparse
import os.path
import sys
import math

def tecplot_to_vtk_camera(psi_angle, theta_angle, viewer_position, view_width, center_of_rotation=None):
    """
    Convert Tecplot 360 camera parameters to VTK camera parameters.
    
    Parameters:
    -----------
    psi_angle : float
        Vertical rotation angle (elevation) in degrees.
    theta_angle : float
        Horizontal rotation angle (azimuth) in degrees.
    viewer_position : tuple or list
        (X, Y, Z) coordinates of the camera position.
    view_width : float
        Width of the viewing window.
    center_of_rotation : tuple or list, optional
        (X, Y, Z) coordinates of the center of rotation (focus point).
        If None, assumes (0, 0, 0).
        
    Returns:
    --------
    dict
        A dictionary containing VTK camera parameters:
        - position: camera position (x, y, z)
        - focal_point: camera focal point (x, y, z)
        - view_up: camera up vector (x, y, z)
        - view_angle: camera view angle in degrees
    """
    # Default center of rotation (focal point) is (0, 0, 0) if not specified
    if center_of_rotation is None:
        center_of_rotation = (0.0, 0.0, 0.0)
    
    # Convert angles from degrees to radians
    psi_rad = math.radians(psi_angle)
    theta_rad = math.radians(theta_angle)
    
    # In Tecplot, the camera position is given directly
    position = viewer_position
    
    # The focal point is the center of rotation
    focal_point = center_of_rotation
    
    # Calculate view up vector based on Tecplot's angle conventions
    # In Tecplot, psi is elevation and theta is azimuth
    # The view up vector points in the direction of the "up" on the screen
    # Using standard spherical to Cartesian conversion with adjusted angles
    view_up_x = -math.sin(theta_rad) * math.sin(psi_rad)
    view_up_y = -math.cos(theta_rad) * math.sin(psi_rad)
    view_up_z = math.cos(psi_rad)
    view_up = (view_up_x, view_up_y, view_up_z)
    
    # Convert view width to view angle
    # In VTK, the view angle is the vertical field of view in degrees
    # Calculate distance from camera to focal point
    dx = position[0] - focal_point[0]
    dy = position[1] - focal_point[1]
    dz = position[2] - focal_point[2]
    distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Convert view width to view angle (approximation)
    # Assuming view_width represents the width of the viewing window in world coordinates
    # VTK uses view angle (field of view) instead
    view_angle = 2 * math.degrees(math.atan2(view_width/2, distance))
    
    return {
        'position': position,
        'focal_point': focal_point,
        'view_up': view_up,
        'view_angle': view_angle
    }

def apply_to_vtk_camera(camera, tecplot_params):
    """
    Apply the converted Tecplot parameters to a VTK camera object.
    
    Parameters:
    -----------
    camera : vtkCamera
        The VTK camera object to modify.
    tecplot_params : dict
        The output of tecplot_to_vtk_camera function.
    """
    camera.SetPosition(*tecplot_params['position'])
    camera.SetFocalPoint(*tecplot_params['focal_point'])
    camera.SetViewUp(*tecplot_params['view_up'])
    camera.SetViewAngle(tecplot_params['view_angle'])

def parse_tecplot_layout(filename):
    """
    Parse a Tecplot layout file (.lay) to extract THREEDVIEW parameters.
    
    Parameters:
    -----------
    filename : str
        Path to the Tecplot layout file.
        
    Returns:
    --------
    dict
        Dictionary containing the parsed parameters.
    """
    params = {
        'psi_angle': None,
        'theta_angle': None,
        'viewer_position': None,
        'view_width': None,
        'center_of_rotation': None
    }
    
    try:
        with open(filename, 'r') as f:
            in_threedview = False
            in_viewer_position = False
            in_center_rotation = False
            
            for line in f:
                line = line.strip()
                
                if line.startswith('$!THREEDVIEW'):
                    in_threedview = True
                    continue
                    
                if not in_threedview:
                    continue
                    
                if line.startswith('PSIANGLE'):
                    parts = line.split('=')
                    if len(parts) > 1:
                        params['psi_angle'] = float(parts[1].strip())
                        
                elif line.startswith('THETAANGLE'):
                    parts = line.split('=')
                    if len(parts) > 1:
                        params['theta_angle'] = float(parts[1].strip())
                        
                elif line.startswith('VIEWERPOSITION'):
                    in_viewer_position = True
                    continue
                    
                elif in_viewer_position and line.startswith('X'):
                    parts = line.split('=')
                    if len(parts) > 1:
                        x = float(parts[1].strip())
                        if params['viewer_position'] is None:
                            params['viewer_position'] = [0, 0, 0]
                        params['viewer_position'][0] = x
                        
                elif in_viewer_position and line.startswith('Y'):
                    parts = line.split('=')
                    if len(parts) > 1:
                        y = float(parts[1].strip())
                        if params['viewer_position'] is None:
                            params['viewer_position'] = [0, 0, 0]
                        params['viewer_position'][1] = y
                        
                elif in_viewer_position and line.startswith('Z'):
                    parts = line.split('=')
                    if len(parts) > 1:
                        z = float(parts[1].strip())
                        if params['viewer_position'] is None:
                            params['viewer_position'] = [0, 0, 0]
                        params['viewer_position'][2] = z
                        
                elif line.startswith('}') and in_viewer_position:
                    in_viewer_position = False
                    
                elif line.startswith('CENTEROFROTATION'):
                    in_center_rotation = True
                    continue
                    
                elif in_center_rotation and line.startswith('X'):
                    parts = line.split('=')
                    if len(parts) > 1:
                        x = float(parts[1].strip())
                        if params['center_of_rotation'] is None:
                            params['center_of_rotation'] = [0, 0, 0]
                        params['center_of_rotation'][0] = x
                        
                elif in_center_rotation and line.startswith('Y'):
                    parts = line.split('=')
                    if len(parts) > 1:
                        y = float(parts[1].strip())
                        if params['center_of_rotation'] is None:
                            params['center_of_rotation'] = [0, 0, 0]
                        params['center_of_rotation'][1] = y
                        
                elif in_center_rotation and line.startswith('Z'):
                    parts = line.split('=')
                    if len(parts) > 1:
                        z = float(parts[1].strip())
                        if params['center_of_rotation'] is None:
                            params['center_of_rotation'] = [0, 0, 0]
                        params['center_of_rotation'][2] = z
                        
                elif line.startswith('}') and in_center_rotation:
                    in_center_rotation = False
                
                elif line.startswith('VIEWWIDTH'):
                    parts = line.split('=')
                    if len(parts) > 1:
                        params['view_width'] = float(parts[1].strip())
                        
                # If there's a blank line or we reach a new section, exit the THREEDVIEW section
                elif line.startswith('$!') and not line.startswith('$!THREEDVIEW'):
                    in_threedview = False
                    
        if params['psi_angle'] is None or params['theta_angle'] is None or params['viewer_position'] is None or params['view_width'] is None:
            print("Warning: Some camera parameters were not found in the layout file.")
            
        return params
        
    except Exception as e:
        print(f"Error parsing Tecplot layout file: {e}")
        return None

def process_tecplot_file(input_file, layout_file):
    """
    Process a Tecplot data file and visualize it using VTK.
    
    Parameters:
    -----------
    input_file : str
        Path to the Tecplot .dat file
    layout_file : str
        Path to the Tecplot .lay layout file with camera parameters
    """
    print(f"Processing data file: {input_file}")
    print(f"Using layout file: {layout_file}")
    
    # Parse layout file
    camera_params = None
    print(f"Reading camera settings from layout file: {layout_file}")
    tecplot_params = parse_tecplot_layout(layout_file)
    if not tecplot_params:
        print("Error: Failed to extract camera parameters from layout file.")
        sys.exit(1)
    
    # Convert Tecplot parameters to VTK camera parameters
    camera_params = tecplot_to_vtk_camera(
        tecplot_params['psi_angle'],
        tecplot_params['theta_angle'],
        tecplot_params['viewer_position'],
        tecplot_params['view_width'],
        tecplot_params['center_of_rotation']
    )
    print("Camera parameters found in layout file:")
    print(f"  Position: {camera_params['position']}")
    print(f"  Focal Point: {camera_params['focal_point']}")
    print(f"  View Up: {camera_params['view_up']}")
    print(f"  View Angle: {camera_params['view_angle']}")
    
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

    # 4. Set up the render window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(800, 800)  # Matching the layout's WIDTH = 8, HEIGHT = 8 but in pixels

    # 5. Set up camera
    camera = renderer.GetActiveCamera()
    
    # Apply camera parameters from layout file
    apply_to_vtk_camera(camera, camera_params)

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
    
    # Add required arguments
    parser.add_argument('input_file', type=str, help='Path to the input Tecplot .dat file')
    parser.add_argument('layout_file', type=str, help='Path to the Tecplot .lay layout file for camera settings')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Make sure both files exist
    if not os.path.isfile(args.input_file):
        print(f"Error: Input data file '{args.input_file}' does not exist.")
        sys.exit(1)
        
    if not os.path.isfile(args.layout_file):
        print(f"Error: Layout file '{args.layout_file}' does not exist.")
        sys.exit(1)
    
    # Process the Tecplot file
    process_tecplot_file(args.input_file, args.layout_file)

if __name__ == '__main__':
    main()