import vtk
import numpy as np
import os.path
import sys
import math

# Import the magnification modules
from . import magnification


class WorkspaceBar:
    """
    A workspace indicator bar in the style of i3wm.
    Displays workspace buttons at the bottom of the VTK window.
    """
    def __init__(self, renderer, num_workspaces=10):
        self.renderer = renderer
        self.num_workspaces = num_workspaces
        self.current_workspace = 0
        self.actors = []
        self.text_actors = []
        self.background_actors = []

        # i3wm color scheme
        self.focused_bg = (0x28/255.0, 0x55/255.0, 0x77/255.0)  # #285577
        self.focused_text = (1.0, 1.0, 1.0)  # #ffffff
        self.inactive_bg = (0x22/255.0, 0x22/255.0, 0x22/255.0)  # #222222
        self.inactive_text = (0x88/255.0, 0x88/255.0, 0x88/255.0)  # #888888
        self.border_color = (0x4c/255.0, 0x78/255.0, 0x99/255.0)  # #4c7899

        self._create_bar()

    def _create_bar(self):
        """Create the workspace bar with buttons for each workspace."""
        button_width = 0.06
        button_height = 0.04
        button_spacing = 0.002
        bar_height = 0.05
        start_x = 0.01
        start_y = 0.01

        for i in range(self.num_workspaces):
            x_pos = start_x + i * (button_width + button_spacing)

            # Create background for button
            bg_actor = vtk.vtkActor2D()
            bg_mapper = vtk.vtkPolyDataMapper2D()

            # Create a rectangle for the button background
            points = vtk.vtkPoints()
            points.InsertNextPoint(x_pos, start_y, 0)
            points.InsertNextPoint(x_pos + button_width, start_y, 0)
            points.InsertNextPoint(x_pos + button_width, start_y + button_height, 0)
            points.InsertNextPoint(x_pos, start_y + button_height, 0)

            quad = vtk.vtkCellArray()
            quad.InsertNextCell(4)
            quad.InsertCellPoint(0)
            quad.InsertCellPoint(1)
            quad.InsertCellPoint(2)
            quad.InsertCellPoint(3)

            poly_data = vtk.vtkPolyData()
            poly_data.SetPoints(points)
            poly_data.SetPolys(quad)

            bg_mapper.SetInputData(poly_data)
            bg_actor.SetMapper(bg_mapper)

            # Set initial color based on whether this is the current workspace
            if i == self.current_workspace:
                bg_actor.GetProperty().SetColor(*self.focused_bg)
            else:
                bg_actor.GetProperty().SetColor(*self.inactive_bg)

            self.background_actors.append(bg_actor)
            self.renderer.AddActor2D(bg_actor)

            # Create text label for workspace number
            text_actor = vtk.vtkTextActor()
            text_actor.SetInput(str(i))
            text_actor.GetTextProperty().SetFontSize(18)
            text_actor.GetTextProperty().SetBold(True)
            text_actor.GetTextProperty().SetJustificationToCentered()
            text_actor.GetTextProperty().SetVerticalJustificationToCentered()

            # Position text in center of button
            # Convert normalized coords to display coords for text positioning
            text_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
            text_actor.GetPositionCoordinate().SetValue(
                x_pos + button_width/2,
                start_y + button_height/2
            )

            # Set text color based on workspace state
            if i == self.current_workspace:
                text_actor.GetTextProperty().SetColor(*self.focused_text)
            else:
                text_actor.GetTextProperty().SetColor(*self.inactive_text)

            self.text_actors.append(text_actor)
            self.renderer.AddActor2D(text_actor)

    def switch_workspace(self, workspace_index):
        """Switch to a different workspace and update the visual appearance."""
        if workspace_index < 0 or workspace_index >= self.num_workspaces:
            return

        old_workspace = self.current_workspace
        self.current_workspace = workspace_index

        # Update old workspace button to inactive style
        self.background_actors[old_workspace].GetProperty().SetColor(*self.inactive_bg)
        self.text_actors[old_workspace].GetTextProperty().SetColor(*self.inactive_text)

        # Update new workspace button to focused style
        self.background_actors[workspace_index].GetProperty().SetColor(*self.focused_bg)
        self.text_actors[workspace_index].GetTextProperty().SetColor(*self.focused_text)


def extract_variable_names(multiblock_data):
    """
    Extract variable names from VTK multiblock data.

    Parameters:
    -----------
    multiblock_data : vtkMultiBlockDataSet
        The multiblock dataset from Tecplot reader

    Returns:
    --------
    list
        List of variable names (e.g., ["x", "y", "z", "stress_vm", ...])
    """
    # Get the first block to extract variable names
    for i in range(multiblock_data.GetNumberOfBlocks()):
        block = multiblock_data.GetBlock(i)
        if block is not None:
            point_data = block.GetPointData()
            num_arrays = point_data.GetNumberOfArrays()
            var_names = []
            for j in range(num_arrays):
                var_names.append(point_data.GetArrayName(j))

            # VTK doesn't include x, y, z in point data arrays (they're used as coordinates)
            # But Tecplot VAR indices start at 1 and include x, y, z
            # So we need to prepend x, y, z to match Tecplot VAR indexing:
            # VAR=1 -> x, VAR=2 -> y, VAR=3 -> z, VAR=4 -> first scalar variable, etc.
            full_var_names = ["x", "y", "z"] + var_names
            print(f"Found {len(full_var_names)} variables (including x,y,z coordinates): {full_var_names}")
            return full_var_names

    # Fallback if no blocks found
    print("Warning: No blocks found in data, using default variable names")
    return ["x", "y", "z", "u"]


def create_contour_group_actors(block, contour_group, lut, variable_names, is_active=True):
    """
    Create VTK actors for all contours in a contour group.

    Parameters:
    -----------
    block : vtkDataSet
        The data block to create contours from
    contour_group : dict
        Dictionary with all GLOBALCONTOUR parameters:
        - 'values': contour level values
        - 'var': variable index
        - 'cmin'/'cmax': color map range
        - 'colorcutoff_min'/'colorcutoff_max': color cutoff range
        - etc.
    lut : vtkLookupTable
        Color lookup table for mapping
    variable_names : list
        List of variable names from data file (indexed from 1, where 1=x, 2=y, 3=z, ...)
    is_active : bool
        Whether this is an active block

    Returns:
    --------
    list
        List of actors for all contours in the group plus edge and surface actors
    """
    actors = []
    contour_values = contour_group['values']

    # Use COLORCUTOFF for scalar range if available, otherwise use CMIN/CMAX
    # COLORCUTOFF defines the visible range, while CMIN/CMAX define the full color map range
    if contour_group.get('colorcutoff_min') is not None and contour_group.get('colorcutoff_max') is not None:
        scalar_range = (contour_group['colorcutoff_min'], contour_group['colorcutoff_max'])
    else:
        scalar_range = (contour_group.get('cmin', 0), contour_group.get('cmax', 1))

    # Get variable name from VAR index (1-based: 1=first variable, etc.)
    var_index = contour_group.get('var', 1)
    if var_index and 0 < var_index <= len(variable_names):
        var_name = variable_names[var_index - 1]  # Convert to 0-based index
    else:
        # Fallback to "u" if VAR index is invalid
        var_name = "u"
        print(f"Warning: Invalid VAR index {var_index}, using fallback variable 'u'")

    print(f"Using variable '{var_name}' for VAR={var_index}")

    # Create contour filter with ALL values from the group
    contour = vtk.vtkContourFilter()
    contour.SetInputData(block)
    contour.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, var_name)

    # Add all contour values
    for idx, value in enumerate(contour_values):
        contour.SetValue(idx, value)

    contour.Update()

    # Create a mapper for the contour output
    contour_mapper = vtk.vtkPolyDataMapper()
    contour_mapper.SetInputConnection(contour.GetOutputPort())
    contour_mapper.SetScalarRange(scalar_range[0], scalar_range[1])
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
    edge_actor.GetProperty().SetRepresentationToWireframe()
    edge_actor.GetProperty().SetColor(0, 0, 0)  # Black
    edge_actor.GetProperty().SetLineWidth(0.1)
    actors.append(edge_actor)

    # Create surface mapper for flood contour
    surface_mapper = vtk.vtkDataSetMapper()
    surface_mapper.SetInputData(block)
    surface_mapper.SetScalarRange(scalar_range[0], scalar_range[1])
    surface_mapper.SetScalarModeToUsePointFieldData()
    surface_mapper.SelectColorArray(var_name)
    surface_mapper.SetLookupTable(lut)

    # Create surface actor
    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(surface_mapper)
    surface_actor.GetProperty().SetInterpolationToGouraud()
    actors.append(surface_actor)

    return actors


class KeyPressCallback:
    """
    Custom callback class for handling keyboard events in VTK.
    This allows overriding default VTK key bindings with custom behavior.
    """
    def __init__(self, interactor, workspace_bar, renderer, render_window,
                 blocks_data, contour_groups, lut, active_fieldmaps,
                 inactive_actors, scalar_bar, variable_names):
        self.interactor = interactor
        self.workspace_bar = workspace_bar
        self.renderer = renderer
        self.render_window = render_window
        self.blocks_data = blocks_data  # List of (block, is_active) tuples
        self.contour_groups = contour_groups  # List of contour group dicts
        self.lut = lut
        self.active_fieldmaps = active_fieldmaps
        self.inactive_actors = inactive_actors
        self.scalar_bar = scalar_bar  # Scalar bar for legend
        self.variable_names = variable_names  # List of variable names from data
        self.current_actors = []  # Currently visible actors

    def dispose_current_contour_group(self):
        """Remove and dispose current contour group actors."""
        for actor in self.current_actors:
            self.renderer.RemoveActor(actor)
        self.current_actors.clear()

    def create_and_show_contour_group(self, workspace_num):
        """Create and display all contours for the specified workspace (contour group)."""
        contour_group = self.contour_groups[workspace_num]

        print(f"Creating CONTOURGROUP {workspace_num} with {len(contour_group['values'])} contours")

        # Update scalar bar for this contour group
        # Update the lookup table range
        if contour_group.get('colorcutoff_min') is not None and contour_group.get('colorcutoff_max') is not None:
            scalar_range = (contour_group['colorcutoff_min'], contour_group['colorcutoff_max'])
        else:
            scalar_range = (contour_group.get('cmin', 0), contour_group.get('cmax', 1))

        self.lut.SetTableRange(scalar_range[0], scalar_range[1])

        # Update scalar bar title (clear if no header defined)
        header = contour_group.get('header')
        if header:
            self.scalar_bar.SetTitle(header)
        else:
            self.scalar_bar.SetTitle("")  # Clear title if no header

        # Create actors for each block
        for block, is_active in self.blocks_data:
            if block is None:
                continue

            if is_active:
                actors = create_contour_group_actors(block, contour_group, self.lut, self.variable_names, is_active)
                for actor in actors:
                    self.renderer.AddActor(actor)
                    self.current_actors.append(actor)

    def __call__(self, obj, event):
        # Get the key that was pressed
        key = self.interactor.GetKeySym()

        # Handle numeric keys 0-9
        if key in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
            workspace_num = int(key)

            # Check if this workspace exists
            if workspace_num < len(self.contour_groups):
                print(f"Switching to workspace {workspace_num}")

                # Dispose old contour group
                self.dispose_current_contour_group()

                # Create and show new contour group
                self.create_and_show_contour_group(workspace_num)

                # Update workspace bar
                self.workspace_bar.switch_workspace(workspace_num)

                # Render the changes
                self.render_window.Render()
            else:
                print(f"Workspace {workspace_num} does not exist (only {len(self.contour_groups)} workspaces available)")

            # Prevent the default VTK behavior by not calling OnKeyPress
            return

        # For all other keys, let the default interactor style handle them
        # Get the current interactor style and call its OnKeyPress method
        style = self.interactor.GetInteractorStyle()
        if style and hasattr(style, 'OnKeyPress'):
            style.OnKeyPress()
            style.OnChar()  # Also call OnChar to ensure proper handling


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
    Parse a Tecplot layout file (.lay) to extract THREEDVIEW parameters, data file path, and CONTOURGROUPS.

    Parameters:
    -----------
    filename : str
        Path to the Tecplot layout file.

    Returns:
    --------
    dict
        Dictionary containing the parsed parameters including contour_groups list.
    """
    params = {
        'psi_angle': None,
        'theta_angle': None,
        'viewer_position': None,
        'view_width': None,
        'active_fieldmaps': None,
        'data_file': None,
        'legend_header': None,
        'contour_groups': [],  # List of contour groups, each with var, values, cmin, cmax
    }
    
    try:
        with open(filename, 'r') as f:
            in_threedview = False
            in_viewer_position = False
            in_center_rotation = False
            in_rotate_origin = False
            in_globalcontour = False
            in_legend = False
            in_contourlevels = False
            in_rawdata = False

            # Store all GLOBALCONTOUR sections by group number
            globalcontour_dict = {}
            current_globalcontour = {}
            current_contourlevels = {}
            rawdata_lines = []
            rawdata_count = 0
            rawdata_expected = 0

            # Additional states for nested structures
            in_colorcutoff = False
            in_colormapfilter = False
            in_continuouscolor = False
            in_labels = False

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

                # Parse GLOBALCONTOUR sections
                if line.startswith('$!GLOBALCONTOUR'):
                    # Save the previous GLOBALCONTOUR if it exists and is a different group
                    if in_globalcontour and current_globalcontour.get('group') is not None:
                        prev_group = current_globalcontour['group']
                        globalcontour_dict[prev_group] = current_globalcontour.copy()

                    # Extract contour group number
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            group_num = int(parts[1])

                            # If this group already exists, merge with existing data
                            # Otherwise create new entry
                            if group_num in globalcontour_dict:
                                # Continue with existing data to merge properties
                                current_globalcontour = globalcontour_dict[group_num].copy()
                            else:
                                # Create new entry with all None values
                                current_globalcontour = {
                                    'group': group_num,
                                    'var': None,
                                    'cmin': None,
                                    'cmax': None,
                                    'header': None,
                                    'defnumlevels': None,
                                    'autolevelskip': None,
                                    'legend_show': None,
                                    'legend_xypos': None,
                                    'colorcutoff_min': None,
                                    'colorcutoff_max': None
                                }

                            in_globalcontour = True
                            # Reset nested states
                            in_legend = False
                            in_labels = False
                            in_colorcutoff = False
                            in_colormapfilter = False
                            in_continuouscolor = False
                            # Only save legend header from first contour group
                            if group_num == 1 and params['legend_header'] is None:
                                # Will be parsed in the section below
                                pass
                        except:
                            pass
                    continue

                # Process GLOBALCONTOUR section
                if in_globalcontour:
                    if line.startswith('VAR'):
                        parts = line.split('=')
                        if len(parts) > 1:
                            current_globalcontour['var'] = int(parts[1].strip())
                    elif line.startswith('DEFNUMLEVELS'):
                        parts = line.split('=')
                        if len(parts) > 1:
                            current_globalcontour['defnumlevels'] = int(parts[1].strip())
                    elif line.startswith('LABELS'):
                        in_labels = True
                        continue
                    elif in_labels:
                        if line.startswith('AUTOLEVELSKIP'):
                            parts = line.split('=')
                            if len(parts) > 1:
                                current_globalcontour['autolevelskip'] = int(parts[1].strip())
                        elif line.startswith('}'):
                            in_labels = False
                    elif line.startswith('LEGEND'):
                        in_legend = True
                        continue
                    elif in_legend:
                        if line.startswith('SHOW'):
                            parts = line.split('=')
                            if len(parts) > 1:
                                current_globalcontour['legend_show'] = parts[1].strip().upper() == 'YES'
                        elif line.startswith('HEADER'):
                            parts = line.split('=')
                            if len(parts) > 1:
                                header = parts[1].strip().strip('"')
                                current_globalcontour['header'] = header
                                # Save first contour group header as legend_header
                                if current_globalcontour.get('group') == 1 and params['legend_header'] is None:
                                    params['legend_header'] = header
                                    print(f"Found legend header: {params['legend_header']}")
                        elif line.startswith('X') and '=' in line:
                            parts = line.split('=')
                            if len(parts) > 1 and current_globalcontour['legend_xypos'] is None:
                                current_globalcontour['legend_xypos'] = [float(parts[1].strip()), 0]
                        elif line.startswith('Y') and '=' in line:
                            parts = line.split('=')
                            if len(parts) > 1 and current_globalcontour['legend_xypos']:
                                current_globalcontour['legend_xypos'][1] = float(parts[1].strip())
                        elif line.startswith('}'):
                            in_legend = False
                    elif line.startswith('COLORCUTOFF'):
                        in_colorcutoff = True
                        continue
                    elif in_colorcutoff:
                        if line.startswith('RANGEMIN'):
                            parts = line.split('=')
                            if len(parts) > 1:
                                current_globalcontour['colorcutoff_min'] = float(parts[1].strip())
                        elif line.startswith('RANGEMAX'):
                            parts = line.split('=')
                            if len(parts) > 1:
                                current_globalcontour['colorcutoff_max'] = float(parts[1].strip())
                        elif line.startswith('}'):
                            in_colorcutoff = False
                    elif line.startswith('COLORMAPFILTER'):
                        in_colormapfilter = True
                        continue
                    elif in_colormapfilter:
                        if line.startswith('CONTINUOUSCOLOR'):
                            in_continuouscolor = True
                            continue
                        elif in_continuouscolor:
                            if line.startswith('CMIN'):
                                parts = line.split('=')
                                if len(parts) > 1:
                                    current_globalcontour['cmin'] = float(parts[1].strip())
                            elif line.startswith('CMAX'):
                                parts = line.split('=')
                                if len(parts) > 1:
                                    current_globalcontour['cmax'] = float(parts[1].strip())
                            elif line.startswith('}'):
                                in_continuouscolor = False
                        elif line.startswith('}') and not in_continuouscolor:
                            in_colormapfilter = False
                    elif line.startswith('$!') and not line.startswith('$!GLOBALCONTOUR'):
                        # End of GLOBALCONTOUR section - store it in the dict
                        if current_globalcontour.get('group') is not None:
                            globalcontour_dict[current_globalcontour['group']] = current_globalcontour.copy()
                        in_globalcontour = False
                        in_legend = False
                        in_labels = False
                        in_colorcutoff = False
                        in_colormapfilter = False
                        in_continuouscolor = False

                # Parse CONTOURLEVELS sections
                if line.startswith('$!CONTOURLEVELS'):
                    in_contourlevels = True
                    current_contourlevels = {}
                    continue

                if in_contourlevels:
                    if line.startswith('CONTOURGROUP'):
                        parts = line.split('=')
                        if len(parts) > 1:
                            current_contourlevels['group'] = int(parts[1].strip())
                    elif line.startswith('RAWDATA'):
                        in_rawdata = True
                        rawdata_lines = []
                        rawdata_count = 0
                        continue
                    elif in_rawdata:
                        # Try to parse as number
                        try:
                            if rawdata_count == 0:
                                # First line is the count
                                rawdata_expected = int(line)
                                rawdata_count += 1
                            else:
                                # Subsequent lines are values
                                rawdata_lines.append(float(line))
                                if len(rawdata_lines) >= rawdata_expected:
                                    # Finished reading rawdata
                                    current_contourlevels['values'] = rawdata_lines
                                    in_rawdata = False
                                    in_contourlevels = False

                                    # Combine with GLOBALCONTOUR info
                                    group_num = current_contourlevels.get('group')
                                    if group_num:
                                        # Find matching GLOBALCONTOUR data
                                        globalcontour_data = globalcontour_dict.get(group_num, {})
                                        group_data = {
                                            'group': group_num,
                                            'values': current_contourlevels['values'],
                                            'var': globalcontour_data.get('var'),
                                            'cmin': globalcontour_data.get('cmin'),
                                            'cmax': globalcontour_data.get('cmax'),
                                            'header': globalcontour_data.get('header'),
                                            'defnumlevels': globalcontour_data.get('defnumlevels'),
                                            'autolevelskip': globalcontour_data.get('autolevelskip'),
                                            'legend_show': globalcontour_data.get('legend_show'),
                                            'legend_xypos': globalcontour_data.get('legend_xypos'),
                                            'colorcutoff_min': globalcontour_data.get('colorcutoff_min'),
                                            'colorcutoff_max': globalcontour_data.get('colorcutoff_max')
                                        }
                                        params['contour_groups'].append(group_data)
                                        print(f"Found CONTOURGROUP {group_num} with {len(rawdata_lines)} levels")
                        except ValueError:
                            # Not a number, might be end of section
                            if line.startswith('$!'):
                                in_rawdata = False
                                in_contourlevels = False
                
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
                            
        # After parsing, check for GLOBALCONTOUR sections without CONTOURLEVELS
        # and generate contour values automatically
        existing_groups = {g['group'] for g in params['contour_groups']}

        for group_num, globalcontour_data in globalcontour_dict.items():
            if group_num not in existing_groups:
                # This GLOBALCONTOUR has no corresponding CONTOURLEVELS
                # Generate contour values automatically
                cmin = globalcontour_data.get('cmin')
                cmax = globalcontour_data.get('cmax')

                # Skip groups without CMIN/CMAX
                if cmin is None or cmax is None:
                    print(f"Warning: CONTOURGROUP {group_num} has no CMIN/CMAX, skipping")
                    continue

                defnumlevels = globalcontour_data.get('defnumlevels', 10)

                if defnumlevels and defnumlevels > 0:
                    # Generate equally spaced contour values
                    values = []
                    for i in range(defnumlevels):
                        value = cmin + (cmax - cmin) * i / (defnumlevels - 1) if defnumlevels > 1 else cmin
                        values.append(value)
                else:
                    # No DEFNUMLEVELS specified, use default 10 levels
                    values = [cmin + (cmax - cmin) * i / 9 for i in range(10)]

                group_data = {
                    'group': group_num,
                    'values': values,
                    'var': globalcontour_data.get('var'),
                    'cmin': cmin,
                    'cmax': cmax,
                    'header': globalcontour_data.get('header'),
                    'defnumlevels': defnumlevels,
                    'autolevelskip': globalcontour_data.get('autolevelskip'),
                    'legend_show': globalcontour_data.get('legend_show'),
                    'legend_xypos': globalcontour_data.get('legend_xypos'),
                    'colorcutoff_min': globalcontour_data.get('colorcutoff_min'),
                    'colorcutoff_max': globalcontour_data.get('colorcutoff_max')
                }
                params['contour_groups'].append(group_data)
                print(f"Found CONTOURGROUP {group_num} with {len(values)} auto-generated levels")

        # Sort contour groups by group number
        params['contour_groups'].sort(key=lambda g: g['group'])

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

    # Extract variable names from the data
    variable_names = extract_variable_names(multiblock_data)

    # Create a lookup table for color mapping (will be used for both active and inactive blocks)
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfColors(256)
    lut.SetHueRange(0.667, 0.0)  # Reverse the hue range (red for high values, blue for low)
    lut.SetSaturationRange(1.0, 1.0)
    lut.SetValueRange(1.0, 1.0)
    lut.Build()

    # Get contour groups from parsed layout parameters
    contour_groups = tecplot_params.get('contour_groups', [])

    if not contour_groups:
        print("Warning: No CONTOURGROUP data found in layout file, using default values")
        # Fallback to default values if no contour groups found
        contour_groups = [{
            'group': 1,
            'values': [
                0.335123736549, 0.344903259539, 0.35468278253, 0.364462305521,
                0.374241828511, 0.384021351502, 0.393800874493, 0.403580397483,
                0.413359920474, 0.423139443464, 0.432918966455, 0.442698489446,
                0.452478012436, 0.462257535427, 0.472037058418, 0.481816581408,
                0.491596104399, 0.50137562739, 0.51115515038
            ],
            'var': 5,
            'cmin': 0.325344213558,
            'cmax': 0.520934673371,
            'header': 'u velocity'
        }]

    # Store block information for lazy contour group creation
    blocks_data = []

    print(f"Preparing {len(contour_groups)} contour group workspaces...")
    for i, group in enumerate(contour_groups):
        print(f"  Workspace {i}: CONTOURGROUP {group['group']} with {len(group['values'])} contours")

    # 2. Process each block in the multi-block dataset and store block info
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

        # Store block data for later use
        blocks_data.append((block, is_active))

    # Create actors ONLY for the first workspace (workspace 0 = first contour group)
    if contour_groups:
        first_group = contour_groups[0]
        print(f"\nCreating initial CONTOURGROUP 0 with {len(first_group['values'])} contours...")
        for block, is_active in blocks_data:
            if block is None:
                continue

            if is_active:
                actors = create_contour_group_actors(block, first_group, lut, variable_names, is_active)
                active_actors.extend(actors)
            else:
                # For inactive blocks, still create surface actors (no contour)
                scalar_range = (first_group.get('cmin', 0), first_group.get('cmax', 1))

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
                surface_mapper.SetScalarRange(scalar_range[0], scalar_range[1])
                surface_mapper.SetScalarModeToUsePointFieldData()
                surface_mapper.SelectColorArray("u")
                surface_mapper.SetLookupTable(lut)

                # Create surface actor
                surface_actor = vtk.vtkActor()
                surface_actor.SetMapper(surface_mapper)
                surface_actor.GetProperty().SetInterpolationToGouraud()

                inactive_actors.extend([edge_actor, surface_actor])

    print(f"Created {len(active_actors)} active actors for initial workspace")

    # 3. Set up the renderer and render window
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(800, 800)

    # 4. Add only the ACTIVE actors from first workspace to renderer
    print(f"Adding {len(active_actors)} active actors to renderer")
    for actor in active_actors:
        renderer.AddActor(actor)

    # 5. Set up the camera
    camera = renderer.GetActiveCamera()
    apply_to_vtk_camera(camera, camera_params)

    # 6. Reset the camera to fit all ACTIVE actors from first workspace
    if active_actors:
        print("Resetting camera to fit active actors")
        renderer.ResetCamera()

        # Apply the selected magnification method if not 'none'
        if magnification_method != 'none':
            print(f"Applying {magnification_method} magnification...")
            if magnification_method == 'statistical':
                magnification.statistical.apply_magnification(renderer, active_actors)
            elif magnification_method == 'variability':
                magnification.variability.apply_magnification(renderer, active_actors)
            elif magnification_method == 'field_range':
                magnification.field_range.apply_magnification(renderer, active_actors)
            elif magnification_method == 'bounding_box':
                magnification.bounding_box.apply_magnification(renderer, active_actors)
            elif magnification_method == 'mesh_refinement':
                magnification.mesh_refinement.apply_magnification(renderer, active_actors)
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
    # Set title only if legend header is defined in layout file
    legend_title = tecplot_params.get('legend_header')
    if legend_title:
        scalar_bar.SetTitle(legend_title)
    scalar_bar.SetNumberOfLabels(10)
    scalar_bar.SetPosition(0.8, 0.1)
    scalar_bar.SetPosition2(0.15, 0.8)
    renderer.AddActor2D(scalar_bar)

    # Create workspace bar in i3wm style
    num_workspaces = min(len(contour_groups), 10)  # Limit to 10 for keys 0-9
    workspace_bar = WorkspaceBar(renderer, num_workspaces)
    print(f"\nWorkspace bar created with {num_workspaces} CONTOURGROUP workspaces")

    # 9. Initialize the interactor
    renderWindowInteractor.Initialize()

    # Remove default CharEvent observers from interactor
    # The interactor processes CharEvent AFTER KeyPressEvent, and this is where
    # the default key handlers (like stereo toggle on '3') are executed.
    # By removing CharEvent observers, we prevent the default handling.
    # Our custom callback will then explicitly call OnChar for non-numeric keys.
    renderWindowInteractor.RemoveObservers('CharEvent')

    # Add custom key press callback to override default VTK key bindings
    key_callback = KeyPressCallback(
        renderWindowInteractor,
        workspace_bar,
        renderer,
        renderWindow,
        blocks_data,
        contour_groups,
        lut,
        tecplot_params.get('active_fieldmaps'),
        inactive_actors,
        scalar_bar,
        variable_names
    )
    # Store initial actors in the callback
    key_callback.current_actors = active_actors.copy()

    renderWindowInteractor.AddObserver('KeyPressEvent', key_callback)
    print("Custom keyboard handler installed:")
    print(f"  - Numeric keys 0-{num_workspaces-1} switch between CONTOURGROUP workspaces")
    print("  - Key '3' stereo mode is disabled")
    print("  - All other keys work as normal (r=reset, w=wireframe, s=surface, etc.)\n")

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
