"""
Bounding box-based view magnification for Tecplot-VTK visualization.

This module implements a magnification method that adjusts the view
based on the bounding box of active blocks with optional padding.
"""

import math
import vtk

def apply_magnification(renderer, active_actors, padding_factor=0.05):
    """
    Apply view magnification using adaptive bounding box approach.
    Adjusts the view based on the bounding box of active blocks.
    
    Parameters:
    -----------
    renderer : vtkRenderer
        The VTK renderer to apply magnification to
    active_actors : list
        List of vtkActor objects representing active blocks
    padding_factor : float, optional
        Factor to add padding around the bounding box (default: 0.05 or 5%)
    """
    # Compute the bounding box of all active actors
    bounds = [float('inf'), float('-inf'), 
              float('inf'), float('-inf'), 
              float('inf'), float('-inf')]  # [xmin, xmax, ymin, ymax, zmin, zmax]
    
    has_valid_bounds = False
    
    for actor in active_actors:
        if not hasattr(actor, 'GetBounds'):
            continue
            
        actor_bounds = actor.GetBounds()
        if not actor_bounds or len(actor_bounds) != 6:
            continue
            
        # Ensure we have valid values before using them
        if all(not math.isinf(b) for b in actor_bounds):
            bounds[0] = min(bounds[0], actor_bounds[0])  # xmin
            bounds[1] = max(bounds[1], actor_bounds[1])  # xmax
            bounds[2] = min(bounds[2], actor_bounds[2])  # ymin
            bounds[3] = max(bounds[3], actor_bounds[3])  # ymax
            bounds[4] = min(bounds[4], actor_bounds[4])  # zmin
            bounds[5] = max(bounds[5], actor_bounds[5])  # zmax
            has_valid_bounds = True
    
    if not has_valid_bounds:
        print("No valid bounding boxes found for bounding-box magnification")
        return
    
    # Calculate the center and dimensions of the bounding box
    center = [
        (bounds[0] + bounds[1]) / 2,
        (bounds[2] + bounds[3]) / 2,
        (bounds[4] + bounds[5]) / 2
    ]
    
    width = bounds[1] - bounds[0]
    height = bounds[3] - bounds[2]
    depth = bounds[5] - bounds[4]
    
    max_dimension = max(width, height, depth)
    
    # Apply padding
    padded_max_dimension = max_dimension * (1 + padding_factor)
    
    # Get the camera and current settings
    camera = renderer.GetActiveCamera()
    position = camera.GetPosition()
    focal_point = camera.GetFocalPoint()
    
    # Calculate direction vector from position to focal point
    direction = [
        focal_point[0] - position[0],
        focal_point[1] - position[1],
        focal_point[2] - position[2]
    ]
    
    # Normalize direction vector
    direction_mag = math.sqrt(sum(d*d for d in direction))
    if direction_mag > 0:
        direction = [d/direction_mag for d in direction]
    
    # Set focal point to the center of the bounding box
    camera.SetFocalPoint(*center)
    
    # Calculate current distance
    current_distance = math.sqrt(sum((p-f)**2 for p, f in zip(position, focal_point)))
    
    # Keep the same distance but change direction to point at the new center
    new_position = [
        center[0] - direction[0] * current_distance,
        center[1] - direction[1] * current_distance, 
        center[2] - direction[2] * current_distance
    ]
    camera.SetPosition(*new_position)
    
    # Calculate zoom factor based on field of view needed to see the entire padded bounding box
    # This is an approximation based on the camera's view angle
    view_angle_rad = math.radians(camera.GetViewAngle())
    ideal_distance = (padded_max_dimension / 2) / math.tan(view_angle_rad / 2)
    zoom_factor = current_distance / ideal_distance if ideal_distance > 0 else 1.0
    # zoom_factor = min(max(zoom_factor, 0.5), 2.0) # Limit to reasonable range
    
    # Apply zoom
    camera.Dolly(zoom_factor)
    
    # Reset clipping range
    renderer.ResetCameraClippingRange()
    
    print(f"Applied bounding box magnification with padding factor: {padding_factor}")
    print(f"Focus center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    print(f"Zoom factor: {zoom_factor:.2f}")
