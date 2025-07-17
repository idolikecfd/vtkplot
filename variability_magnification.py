"""
Variability-based view magnification for Tecplot-VTK visualization.

This module implements a magnification method that focuses on regions
with high data variability (gradients).
"""

import numpy as np
import vtk
import math

def apply_magnification(renderer, active_actors, scalar_name="u", gradient_percentile=0.9):
    """
    Apply view magnification focusing on high data variability regions.
    Identifies areas with high gradient values and zooms to them.
    
    Parameters:
    -----------
    renderer : vtkRenderer
        The VTK renderer to apply magnification to
    active_actors : list
        List of vtkActor objects representing active blocks
    scalar_name : str, optional
        Name of the scalar field to analyze gradients (default: "u")
    gradient_percentile : float, optional
        Percentile threshold for high gradients (default: 0.9 means top 10%)
    """
    # Collect all points with gradient information
    all_points = []
    all_gradient_magnitudes = []
    
    for actor in active_actors:
        if not hasattr(actor, 'GetMapper') or not actor.GetMapper() or not actor.GetMapper().GetInput():
            continue
            
        data = actor.GetMapper().GetInput()
        
        # Skip if no point data
        if not hasattr(data, 'GetPointData') or not data.GetPointData() or not data.GetPointData().GetArray(scalar_name):
            continue
        
        # Calculate gradients
        gradients = vtk.vtkGradientFilter()
        gradients.SetInputData(data)
        gradients.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, scalar_name)
        gradients.SetResultArrayName("GradientMagnitude")
        gradients.ComputeGradientOff()  # Don't compute gradient vector
        gradients.ComputeDivergenceOff()
        gradients.ComputeVorticityOff() 
        gradients.ComputeQCriterionOff()
        gradients.Update()
        
        gradient_output = gradients.GetOutput()
        if not gradient_output:
            continue
            
        gradient_array = gradient_output.GetPointData().GetArray("GradientMagnitude")
        points = gradient_output.GetPoints()
        
        if not gradient_array or not points:
            continue
        
        # Collect gradient magnitudes and points
        for i in range(gradient_array.GetNumberOfTuples()):
            gradient_mag = abs(gradient_array.GetTuple1(i))  # Use absolute value
            all_gradient_magnitudes.append(gradient_mag)
            all_points.append(points.GetPoint(i))
    
    if not all_points or not all_gradient_magnitudes:
        print("No valid gradient data found for variability magnification")
        return
    
    # Find high gradient threshold
    threshold = np.percentile(all_gradient_magnitudes, gradient_percentile * 100)
    
    # Filter points with gradient magnitude above threshold
    filtered_points = [
        all_points[i] 
        for i in range(len(all_points)) 
        if all_gradient_magnitudes[i] >= threshold
    ]
    
    if not filtered_points:
        print("No points above gradient threshold")
        return
    
    # Calculate bounds of high variability region
    x_vals = [p[0] for p in filtered_points]
    y_vals = [p[1] for p in filtered_points]
    z_vals = [p[2] for p in filtered_points]
    
    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    z_min, z_max = min(z_vals), max(z_vals)
    
    # Calculate center of high variability region
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
    
    # Calculate full vs filtered region size for zoom factor
    all_x = [p[0] for p in all_points]
    all_y = [p[1] for p in all_points]
    all_z = [p[2] for p in all_points]
    
    full_size = max(
        max(all_x) - min(all_x),
        max(all_y) - min(all_y),
        max(all_z) - min(all_z)
    )
    
    filtered_size = max(
        x_max - x_min,
        y_max - y_min,
        z_max - z_min
    )
    
    # Calculate appropriate zoom factor
    zoom_factor = 1.2 * (full_size / filtered_size) if filtered_size > 0 else 1.2
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
    
    print(f"Applied variability-based magnification with zoom factor: {zoom_factor:.2f}")
    print(f"Focus center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    print(f"High gradient threshold: {threshold:.6f}")
