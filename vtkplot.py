import vtk
import numpy as np
import argparse
import os.path
import sys
import math

# Import the magnification modules
import statistical_magnification
import variability_magnification
import field_range_magnification
import bounding_box_magnification
import mesh_refinement_magnification

def parse_active_fieldmaps(fieldmaps_str):
    """
    Parse a Tecplot active field maps string into a list of block indices.
    
    Parameters:
    -----------
    fieldmaps_str : str
        String containing the active field maps specification, e.g., "[1,3-5,7]"
        
    Returns:
    --------
    list
        List of integer indices of active field maps (1-based)
    """
    # Remove brackets and spaces
    cleaned_str = fieldmaps_str.strip('[]').replace(' ', '')
    
    # Split by commas
    parts = cleaned_str.split(',')
    
    active_maps = []
    for part in parts:
        if '-' in part:
            # This is a range
            start, end = map(int, part.split('-'))
            active_maps.extend(range(start, end + 1))
        else:
            # This is a single number
            if part:  # Skip empty strings
                active_maps.append(int(part))
    
    return active_maps

def calculate_data_center(input_file, active_fieldmaps=None):
    """
    Calculate the center of the data from a Tecplot .dat file.
    
    Parameters:
    -----------
    input_file : str
        Path to the Tecplot .dat file
    active_fieldmaps : list, optional
        List of active field maps (1-based indices). If provided, only these blocks
        will be considered for the center calculation.
        
    Returns:
    --------
    tuple
        (X, Y, Z) coordinates of the center of the data
    """
    # Create a VTK Tecplot reader
    reader = vtk.vtkTecplotReader()
    reader.SetFileName(input_file)
    reader.Update()
    
    # Get the bounds of the entire dataset
    multiblock_data = reader.GetOutput()
    bounds = [0, 0, 0, 0, 0, 0]  # Initialize [xmin, xmax, ymin, ymax, zmin, zmax]
    
    # Process each block to find overall bounds
    is_first_active_block = True
    
    for i in range(multiblock_data.GetNumberOfBlocks()):
        # Skip this block if it's not in the active field maps
        if active_fieldmaps is not None and (i + 1) not in active_fieldmaps:
            continue
            
        block = multiblock_data.GetBlock(i)
        if block is None:
            continue
        
        block_bounds = list(block.GetBounds())  # Convert tuple to list
        # Update min/max for each dimension
        if is_first_active_block:
            bounds = block_bounds.copy()  # Make a copy of the list
            is_first_active_block = False
        else:
            bounds[0] = min(bounds[0], block_bounds[0])  # xmin
            bounds[1] = max(bounds[1], block_bounds[1])  # xmax
            bounds[2] = min(bounds[2], block_bounds[2])  # ymin
            bounds[3] = max(bounds[3], block_bounds[3])  # ymax
            bounds[4] = min(bounds[4], block_bounds[4])  # zmin
            bounds[5] = max(bounds[5], block_bounds[5])  # zmax
    
    # Calculate center coordinates
    center_x = (bounds[0] + bounds[1]) / 2.0
    center_y = (bounds[2] + bounds[3]) / 2.0
    center_z = (bounds[4] + bounds[5]) / 2.0
    
    return (center_x, center_y, center_z)

def tecplot_to_vtk_camera(psi_angle, theta_angle, viewer_position, view_width, input_file, active_fieldmaps=None):
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
    input_file : str
        Path to the Tecplot .dat file, used to calculate center of
        rotation.
    active_fieldmaps : list, optional
        List of active field maps (1-based indices). If provided, only these blocks
        will be considered for the center calculation.
        
    Returns:
    --------
    dict
        A dictionary containing VTK camera parameters:
        - position: camera position (x, y, z)
        - focal_point: camera focal point (x, y, z)
        - view_up: camera up vector (x, y, z)
        - view_angle: camera view angle in degrees
    """
    # Determine center of rotation (focal point)
    if input_file is not None:
        print("Calculating center from data...")
        if active_fieldmaps:
            print(f"Using active field maps: {active_fieldmaps}")
        center_of_rotation = calculate_data_center(input_file, active_fieldmaps)
        print(f"Calculated center of rotation: {center_of_rotation}")
    else:
        center_of_rotation = (0.0, 0.0, 0.0)
        print("Using default center of rotation: (0, 0, 0)")
    
    # Convert angles from degrees to radians
    psi_rad = math.radians(psi_angle)
    theta_rad = math.radians(theta_angle)
    
    # In Tecplot, the camera position is given directly
    position = viewer_position
    
    # The focal point is the center of rotation
    focal_point = center_of_rotation
    
    # Calculate camera direction vector (from position to focal point)
    direction = [
        focal_point[0] - position[0],
        focal_point[1] - position[1],
        focal_point[2] - position[2]
    ]
    direction_magnitude = math.sqrt(direction[0]**2 + direction[1]**2 + direction[2]**2)
    
    # Normalize direction vector
    if direction_magnitude > 0:
        direction = [
            direction[0]/direction_magnitude,
            direction[1]/direction_magnitude,
            direction[2]/direction_magnitude
        ]
    
    # Calculate view up vector based on Tecplot's angle conventions
    # In Tecplot, psi is elevation and theta is azimuth
    # The view up vector points in the direction of the "up" on the screen
    # Using standard spherical to Cartesian conversion with adjusted angles
    view_up_x = -math.sin(theta_rad) * math.sin(psi_rad)
    view_up_y = -math.cos(theta_rad) * math.sin(psi_rad)
    view_up_z = math.cos(psi_rad)
    view_up = [view_up_x, view_up_y, view_up_z]
    
    # Check if view_up and direction are too parallel
    # Calculate dot product
    dot_product = view_up[0]*direction[0] + view_up[1]*direction[1] + view_up[2]*direction[2]
    
    # If they're nearly parallel (dot product close to 1 or -1), adjust view_up
    if abs(abs(dot_product) - 1.0) < 0.01:
        print("Warning: View-up vector nearly parallel to camera direction. Adjusting...")
        # Choose a different view_up vector
        # If direction is mostly along Z, use Y as view_up
        if abs(direction[2]) > 0.9:
            view_up = [0, 1, 0]
        # Otherwise use Z as view_up
        else:
            view_up = [0, 0, 1]
    
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
    Parse a Tecplot layout file (.lay) to extract THREEDVIEW parameters and data file path.
    
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
        'active_fieldmaps': None,  # New parameter
        'data_file': None,  # Path to the data file
    }
    
    try:
        with open(filename, 'r') as f:
            in_threedview = False
            in_viewer_position = False
            in_center_rotation = False
            in_rotate_origin = False
            
            for line in f:
                line = line.strip()
                
                # Check for data file path
                if line.startswith('$!VarSet') and '|LFDSFN1|' in line:
                    # Extract the file path from the line
                    # Format: $!VarSet |LFDSFN1| = '"filename.dat"'
                    parts = line.split('=', 1)
                    if len(parts) > 1:
                        # Remove quotes and extra spaces
                        data_file_str = parts[1].strip()
                        data_file_str = data_file_str.strip("'\"")
                        params['data_file'] = data_file_str
                        print(f"Found data file: {params['data_file']}")
                
                # Check for ACTIVEFIELDMAPS
                if line.startswith('$!ACTIVEFIELDMAPS'):
                    parts = line.split('=')
                    if len(parts) > 1:
                        fieldmaps_str = parts[1].strip()
                        params['active_fieldmaps'] = parse_active_fieldmaps(fieldmaps_str)
                        print(f"Found active fieldmaps: {params['active_fieldmaps']}")
                
                # Check for THREEDVIEW section
                if line.startswith('$!THREEDVIEW'):
                    in_threedview = True
                    continue
                
                # Process THREEDVIEW section
                if in_threedview:
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
                    
                    elif line.startswith('VIEWWIDTH'):
                        parts = line.split('=')
                        if len(parts) > 1:
                            params['view_width'] = float(parts[1].strip())
                            
                    # If there's a blank line or we reach a new section, exit the THREEDVIEW section
                    elif line.startswith('$!') and not line.startswith('$!THREEDVIEW'):
                        in_threedview = False
                            
        if params['psi_angle'] is None or params['theta_angle'] is None or params['viewer_position'] is None or params['view_width'] is None:
            print("Warning: Some camera parameters were not found in the layout file.")
        
        if params['data_file'] is None:
            print("Warning: Data file path was not found in the layout file.")
            
        return params
        
    except Exception as e:
        print(f"Error parsing Tecplot layout file: {e}")
        return None

def process_tecplot_file(layout_file, magnification_method='mesh_refinement'):
    """
    Process a Tecplot layout file and visualize the referenced data using VTK.
    
    Parameters:
    -----------
    layout_file : str
        Path to the Tecplot .lay layout file with camera parameters and data file reference
    magnification_method : str, optional
        View magnification method to use (default: 'mesh_refinement')
    """
    print(f"Using layout file: {layout_file}")
    print(f"Using magnification method: {magnification_method}")
    
    # Parse layout file
    print(f"Reading settings from layout file: {layout_file}")
    tecplot_params = parse_tecplot_layout(layout_file)
    if not tecplot_params:
        print("Error: Failed to extract parameters from layout file.")
        sys.exit(1)
    
    # Get data file path
    if not tecplot_params['data_file']:
        print("Error: No data file path found in layout file.")
        sys.exit(1)
    
    # Resolve data file path relative to layout file directory
    layout_dir = os.path.dirname(os.path.abspath(layout_file))
    input_file = os.path.join(layout_dir, tecplot_params['data_file'])
    
    if not os.path.exists(input_file):
        # Try as absolute path
        input_file = tecplot_params['data_file']
        if not os.path.exists(input_file):
            print(f"Error: Data file '{tecplot_params['data_file']}' not found.")
            print(f"Tried: {os.path.join(layout_dir, tecplot_params['data_file'])}")
            print(f"And: {tecplot_params['data_file']}")
            sys.exit(1)
    
    print(f"Processing data file: {input_file}")
    
    # Convert Tecplot parameters to VTK camera parameters
    camera_params = tecplot_to_vtk_camera(
        tecplot_params['psi_angle'],
        tecplot_params['theta_angle'],
        tecplot_params['viewer_position'],
        tecplot_params['view_width'],
        input_file,  # Pass the input file for center calculation
        tecplot_params.get('active_fieldmaps')  # Pass active fieldmaps
    )
    print("Camera parameters:")
    print(f"  Position: {camera_params['position']}")
    print(f"  Focal Point: {camera_params['focal_point']}")
    print(f"  View Up: {camera_params['view_up']}")
    print(f"  View Angle: {camera_params['view_angle']}")
    
    # 1. Load the Tecplot data file
    reader = vtk.vtkTecplotReader()
    reader.SetFileName(input_file)
    reader.Update()

    # Initialize actors list for active and inactive blocks
    active_actors = []
    inactive_actors = []

    # Get the multi-block output
    multiblock_data = reader.GetOutput()
    print(f"Data type: {multiblock_data.GetClassName()}")
    print(f"Number of blocks: {multiblock_data.GetNumberOfBlocks()}")

    # Create a lookup table for color mapping (will be used for both active and inactive blocks)
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfColors(256)
    lut.SetHueRange(0.667, 0.0)  # Reverse the hue range (red for high values, blue for low)
    lut.SetSaturationRange(1.0, 1.0)
    lut.SetValueRange(1.0, 1.0)
    lut.Build()
    
    # 2. Process each block in the multi-block dataset
    for i in range(multiblock_data.GetNumberOfBlocks()):
        block = multiblock_data.GetBlock(i)
        if block is None:
            print(f"Block {i+1} is None, skipping")
            continue
        
        # Check if the block is active
        is_active = (tecplot_params.get('active_fieldmaps') is None or 
                    (i + 1) in tecplot_params.get('active_fieldmaps'))
        
        if is_active:
            print(f"Processing active block {i+1}, type: {block.GetClassName()}")
        else:
            print(f"Processing inactive block {i+1}, type: {block.GetClassName()}")
        
        # Apply contour settings from the Tecplot file
        contour = vtk.vtkContourFilter()
        contour.SetInputData(block)
        
        # Set up contour values based on Tecplot levels
        contour_values = [
            0.335123736549, 0.344903259539, 0.35468278253, 0.364462305521,
            0.374241828511, 0.384021351502, 0.393800874493, 0.403580397483,
            0.413359920474, 0.423139443464, 0.432918966455, 0.442698489446,
            0.452478012436, 0.462257535427, 0.472037058418, 0.481816581408,
            0.491596104399, 0.50137562739, 0.51115515038
        ]
        
        # Set contour array
        contour.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "u")
        
        for value in contour_values:
            contour.SetValue(contour_values.index(value), value)
        
        contour.Update()
        
        # Create a mapper for the contour output
        contour_mapper = vtk.vtkPolyDataMapper()
        contour_mapper.SetInputConnection(contour.GetOutputPort())
        contour_mapper.SetScalarRange(0.325344213558, 0.520934673371)  # From CMIN/CMAX in Tecplot layout
        contour_mapper.SetLookupTable(lut)
        
        # Create a contour actor
        contour_actor = vtk.vtkActor()
        contour_actor.SetMapper(contour_mapper)
        
        # Create a mapper for the original dataset (for mesh edges)
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(block)
        
        # Create an actor for edges only
        edge_actor = vtk.vtkActor()
        edge_actor.SetMapper(mapper)
        edge_actor.GetProperty().SetRepresentationToWireframe()
        edge_actor.GetProperty().SetColor(0, 0, 0)  # Black
        edge_actor.GetProperty().SetLineWidth(0.1)
        
        # Create surface mapper for flood contour
        surface_mapper = vtk.vtkDataSetMapper()
        surface_mapper.SetInputData(block)
        surface_mapper.SetScalarRange(0.325344213558, 0.520934673371)
        surface_mapper.SetScalarModeToUsePointFieldData()
        surface_mapper.SelectColorArray("u")
        surface_mapper.SetLookupTable(lut)
        
        # Create surface actor
        surface_actor = vtk.vtkActor()
        surface_actor.SetMapper(surface_mapper)
        surface_actor.GetProperty().SetInterpolationToGouraud()
        
        # Add actors to the appropriate list
        if is_active:
            active_actors.extend([contour_actor, edge_actor, surface_actor])
        else:
            inactive_actors.extend([contour_actor, edge_actor, surface_actor])

    # 3. Set up the renderer and render window
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(800, 800)

    # 4. First add all ACTIVE actors
    print(f"Adding {len(active_actors)} active actors to renderer")
    for actor in active_actors:
        renderer.AddActor(actor)
    
    # 5. Set up the camera
    camera = renderer.GetActiveCamera()
    apply_to_vtk_camera(camera, camera_params)
    
    # 6. Reset the camera to fit all ACTIVE actors
    if active_actors:
        print("Resetting camera to fit active actors")
        renderer.ResetCamera()
        
        # Apply the selected magnification method if not 'none'
        if magnification_method != 'none':
            print(f"Applying {magnification_method} magnification...")
            if magnification_method == 'statistical':
                statistical_magnification.apply_magnification(renderer, active_actors)
            elif magnification_method == 'variability':
                variability_magnification.apply_magnification(renderer, active_actors)
            elif magnification_method == 'field_range':
                field_range_magnification.apply_magnification(renderer, active_actors)
            elif magnification_method == 'bounding_box':
                bounding_box_magnification.apply_magnification(renderer, active_actors)
            elif magnification_method == 'mesh_refinement':
                mesh_refinement_magnification.apply_magnification(renderer, active_actors)
            else:
                raise ValueError(f"Unknown magnification method: {magnification_method}")
    
    # 7. Now add all INACTIVE actors without resetting the camera
    print(f"Adding {len(inactive_actors)} inactive actors to renderer (without resetting camera)")
    for actor in inactive_actors:
        renderer.AddActor(actor)

    # Set background color
    renderer.SetBackground(0.6, 0.3, 0)  # Brown/orange background

    # 8. Set up the interactor
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    style = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(style)

    # Create a color bar / scalar bar for the contour values
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(lut)
    scalar_bar.SetTitle("u velocity")
    scalar_bar.SetNumberOfLabels(10)
    scalar_bar.SetPosition(0.8, 0.1)
    scalar_bar.SetPosition2(0.15, 0.8)
    renderer.AddActor2D(scalar_bar)

    # 9. Initialize the interactor
    renderWindowInteractor.Initialize()

    # Add axis display
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
    # Create the main argument parser
    parser = argparse.ArgumentParser(description='Convert and visualize Tecplot data files using VTK')
    
    # Add required argument
    parser.add_argument('layout_file', type=str, help='Path to the Tecplot .lay layout file')
    
    # Create subparsers for different magnification methods
    subparsers = parser.add_subparsers(dest='magnification', 
                                      title='Magnification methods',
                                      description='Select a view magnification method (default: mesh_refinement):',
                                      help='Method to focus view on important data regions')
    
    # Set default value for magnification if no subcommand is provided
    subparsers.default = 'mesh_refinement'
    
    # Create a parser for the "none" command (no magnification)
    none_parser = subparsers.add_parser('none', help='No magnification, standard VTK camera reset')
    
    # Create a parser for statistical magnification
    statistical_parser = subparsers.add_parser('statistical', 
                                             help='Statistical magnification (exclude outliers)')
    statistical_parser.add_argument('--percentile-cutoff', type=float, default=0.05,
                                  help='Percentile cutoff for outlier removal (default: 0.05)')
    
    # Create a parser for variability magnification
    variability_parser = subparsers.add_parser('variability', 
                                             help='Variability magnification (focus on high gradient regions)')
    variability_parser.add_argument('--scalar-name', type=str, default='u',
                                  help='Scalar field name to analyze gradients (default: u)')
    variability_parser.add_argument('--gradient-percentile', type=float, default=0.9,
                                  help='Percentile threshold for high gradients (default: 0.9)')
    
    # Create a parser for field range magnification
    field_range_parser = subparsers.add_parser('field_range', 
                                             help='Field range magnification (focus on specific value ranges)')
    field_range_parser.add_argument('--scalar-name', type=str, default='u',
                                  help='Scalar field name to analyze (default: u)')
    field_range_parser.add_argument('--min-percentile', type=float, default=0.25,
                                  help='Lower percentile for field range (default: 0.25)')
    field_range_parser.add_argument('--max-percentile', type=float, default=0.75,
                                  help='Upper percentile for field range (default: 0.75)')
    
    # Create a parser for bounding box magnification
    bounding_box_parser = subparsers.add_parser('bounding_box', 
                                              help='Bounding box magnification (with optional padding)')
    bounding_box_parser.add_argument('--padding-factor', type=float, default=0.05,
                                   help='Padding factor around bounding box (default: 0.05)')
    
    # Create a parser for mesh refinement magnification
    mesh_refinement_parser = subparsers.add_parser('mesh_refinement',
                                               help='Mesh refinement magnification (focus on smallest cells) - DEFAULT')
    mesh_refinement_parser.add_argument('--percentile-threshold', type=float, default=0.1,
                                     help='Percentile threshold for cell size (default: 0.1 means smallest 10%)')
    mesh_refinement_parser.add_argument('--min-cells', type=int, default=10,
                                     help='Minimum number of cells for refinement region (default: 10)')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Make sure layout file exists
    if not os.path.isfile(args.layout_file):
        print(f"Error: Layout file '{args.layout_file}' does not exist.")
        sys.exit(1)
    
    # Modify the magnification modules with user-provided parameters
    if args.magnification == 'statistical':
        # Replace the default apply_magnification function with a partial function
        # that includes the user-provided percentile cutoff
        original_apply = statistical_magnification.apply_magnification
        statistical_magnification.apply_magnification = lambda renderer, active_actors: original_apply(
            renderer, active_actors, percentile_cutoff=args.percentile_cutoff)
    
    elif args.magnification == 'variability':
        original_apply = variability_magnification.apply_magnification
        variability_magnification.apply_magnification = lambda renderer, active_actors: original_apply(
            renderer, active_actors, scalar_name=args.scalar_name, gradient_percentile=args.gradient_percentile)
    
    elif args.magnification == 'field_range':
        original_apply = field_range_magnification.apply_magnification
        field_range_magnification.apply_magnification = lambda renderer, active_actors: original_apply(
            renderer, active_actors, scalar_name=args.scalar_name, 
            min_percentile=args.min_percentile, max_percentile=args.max_percentile)
    
    elif args.magnification == 'bounding_box':
        original_apply = bounding_box_magnification.apply_magnification
        bounding_box_magnification.apply_magnification = lambda renderer, active_actors: original_apply(
            renderer, active_actors, padding_factor=args.padding_factor)
            
    elif args.magnification == 'mesh_refinement':
        # For mesh_refinement, check if we have the attributes (from subparser) or use defaults
        original_apply = mesh_refinement_magnification.apply_magnification
        if hasattr(args, 'percentile_threshold'):
            # User specified mesh_refinement subcommand with options
            mesh_refinement_magnification.apply_magnification = lambda renderer, active_actors: original_apply(
                renderer, active_actors, percentile_threshold=args.percentile_threshold, 
                min_cells_required=args.min_cells)
        else:
            # Default case - no subcommand specified, use default parameters
            mesh_refinement_magnification.apply_magnification = lambda renderer, active_actors: original_apply(
                renderer, active_actors, percentile_threshold=0.1, min_cells_required=10)
    
    # Process the Tecplot file with the specified magnification method
    process_tecplot_file(args.layout_file, args.magnification)

if __name__ == '__main__':
    main()