"""
Statistical-based view magnification for Tecplot-VTK visualization.

This module implements a magnification method that focuses on the central
portion of data distribution, excluding statistical outliers.
"""

import numpy as np
import math
import vtk

def apply_magnification(renderer, active_actors, percentile_cutoff=0.05):
    """
    Apply view magnification using statistical approach.
    Focuses on the center 90% of data points, excluding outliers.
    
    Parameters:
    -----------
    renderer : vtkRenderer
        The VTK renderer to apply magnification to
    active_actors : list
        List of vtkActor objects representing active blocks
    percentile_cutoff : float, optional
        Percentile cutoff for outlier removal (default: 0.05 means 5% on each end)
    """
    # Collect points from all active actors
    all_points_x = []
    all_points_y = []
    all_points_z = []
    
    for actor in active_actors:
        if not hasattr(actor, 'GetMapper') or not actor.GetMapper() or not actor.GetMapper().GetInput():
            continue
            
        data = actor.GetMapper().GetInput()
        if not hasattr(data, 'GetPoints'):
            continue
            
        points = data.GetPoints()
        if not points:
            continue
            
        num_points = points.GetNumberOfPoints()
        for i in range(num_points):
            point = points.GetPoint(i)
            all_points_x.append(point[0])
            all_points_y.append(point[1])
            all_points_z.append(point[2])
    
    if not all_points_x:
        print("No valid points found for statistical magnification")
        return
    
    # Calculate percentiles to exclude outliers
    x_min = np.percentile(all_points_x, percentile_cutoff * 100)
    x_max = np.percentile(all_points_x, (1 - percentile_cutoff) * 100)
    y_min = np.percentile(all_points_y, percentile_cutoff * 100)
    y_max = np.percentile(all_points_y, (1 - percentile_cutoff) * 100)
    z_min = np.percentile(all_points_z, percentile_cutoff * 100)
    z_max = np.percentile(all_points_z, (1 - percentile_cutoff) * 100)
    
    # Calculate center
    center = ((x_min + x_max)/2, (y_min + y_max)/2, (z_min + z_max)/2)
    
    # Calculate current vs. adjusted field of view
    camera = renderer.GetActiveCamera()
    position = camera.GetPosition()
    
    # Calculate distance from camera to center
    dist = math.sqrt(sum((p-c)**2 for p, c in zip(position, center)))
    
    # Calculate full range vs. trimmed range
    full_range = max(
        max(all_points_x) - min(all_points_x),
        max(all_points_y) - min(all_points_y),
        max(all_points_z) - min(all_points_z)
    )
    
    trimmed_range = max(
        x_max - x_min,
        y_max - y_min,
        z_max - z_min
    )
    
    # Calculate zoom factor
    zoom_factor = full_range / trimmed_range if trimmed_range > 0 else 1.0
    #zoom_factor = min(max(zoom_factor, 1.0), 2.5) # Limit to reasonable range
    
    print(f"Applying statistical magnification with zoom factor: {zoom_factor:.2f}")
    print(f"Focus center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    
    # Set focal point to the center of the trimmed data
    camera.SetFocalPoint(*center)

    # Apply dolly zoom
    camera.Dolly(zoom_factor)

    renderer.ResetCameraClippingRange()
