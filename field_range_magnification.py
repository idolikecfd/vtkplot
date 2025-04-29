"""
Field value-based view magnification for Tecplot-VTK visualization.

This module implements a magnification method that focuses on regions
where a scalar field falls within a specific value range.
"""

import numpy as np
import math
import vtk

def apply_magnification(renderer, active_actors, scalar_name="u", 
                     min_percentile=0.25, max_percentile=0.75):
    """
    Apply view magnification focusing on specific field value ranges.
    
    Parameters:
    -----------
    renderer : vtkRenderer
        The VTK renderer to apply magnification to
    active_actors : list
        List of vtkActor objects representing active blocks
    scalar_name : str, optional
        Name of the scalar field to analyze (default: "u")
    min_percentile : float, optional
        Lower percentile for field range (default: 0.25)
    max_percentile : float, optional
        Upper percentile for field range (default: 0.75)
    """
    # Collect all field values and their positions
    all_values = []
    all_points = []
    
    for actor in active_actors:
        if not hasattr(actor, 'GetMapper') or not actor.GetMapper() or not actor.GetMapper().GetInput():
            continue
            
        data = actor.GetMapper().GetInput()
        if not hasattr(data, 'GetPoints') or not data.GetPoints():
            continue
            
        points = data.GetPoints()
        
        # Skip if no point data
        if not hasattr(data, 'GetPointData') or not data.GetPointData() or not data.GetPointData().GetArray(scalar_name):
            continue
        
        scalars = data.GetPointData().GetArray(scalar_name)
        
        num_points = points.GetNumberOfPoints()
        for i in range(num_points):
            if i < scalars.GetNumberOfTuples():
                value = scalars.GetTuple1(i)
                point = points.GetPoint(i)
                
                all_values.append(value)
                all_points.append(point)
    
    if not all_values or not all_points:
        print("No valid scalar field data found for field-range magnification")
        return
    
    # Find value range to focus on
    min_threshold = np.percentile(all_values, min_percentile * 100)
    max_threshold = np.percentile(all_values, max_percentile * 100)
    
    # Filter points with values in this range
    filtered_points = [
        all_points[i] 
        for i in range(len(all_points)) 
        if min_threshold <= all_values[i] <= max_threshold
    ]
    
    if not filtered_points:
        print("No points within the specified field value range")
        return
    
    # Calculate bounds of filtered region
    x_vals = [p[0] for p in filtered_points]
    y_vals = [p[1] for p in filtered_points]
    z_vals = [p[2] for p in filtered_points]
    
    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    z_min, z_max = min(z_vals), max(z_vals)
    
    # Calculate center of filtered region
    center = ((x_min + x_max)/2, (y_min + y_max)/2, (z_min + z_max)/2)
    
    # Get current camera settings
    camera = renderer.GetActiveCamera()
    position = camera.GetPosition()
    focal_point = camera.GetFocalPoint()
    
    # Calculate unit direction vector from position to focal point
    direction = [
        focal_point[0] - position[0],
        focal_point[1] - position[1],
        focal_point[2] - position[2]
    ]
    direction_mag = math.sqrt(sum(d*d for d in direction))
    if direction_mag > 0:
        direction = [d/direction_mag for d in direction]
    
    # Calculate current distance
    current_distance = math.sqrt(sum((p-f)**2 for p, f in zip(position, focal_point)))
    
    # Calculate full range vs. filtered range
    all_x = [p[0] for p in all_points]
    all_y = [p[1] for p in all_points]
    all_z = [p[2] for p in all_points]
    
    full_range = max(
        max(all_x) - min(all_x),
        max(all_y) - min(all_y),
        max(all_z) - min(all_z)
    )
    
    filtered_range = max(
        x_max - x_min,
        y_max - y_min,
        z_max - z_min
    )
    
    # Calculate zoom factor - emphasize the region of interest
    zoom_factor = 1.2 * (full_range / filtered_range) if filtered_range > 0 else 1.2
    #zoom_factor = min(max(zoom_factor, 1.0), 2.5)  # Limit to reasonable range
    
    # Set new focal point
    camera.SetFocalPoint(*center)
    
    # Calculate new position based on original direction and distance
    new_position = [
        center[0] - direction[0] * current_distance,
        center[1] - direction[1] * current_distance,
        center[2] - direction[2] * current_distance
    ]
    camera.SetPosition(*new_position)
    
    # Apply zoom
    camera.Dolly(zoom_factor)
    
    # Reset clipping range
    renderer.ResetCameraClippingRange()
    
    print(f"Applied field-range magnification with zoom factor: {zoom_factor:.2f}")
    print(f"Focus center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    print(f"Field value range: {min_threshold:.6f} - {max_threshold:.6f}")
