"""
Mesh refinement-based view magnification for Tecplot-VTK visualization.

This module implements a magnification method that focuses on regions
with the highest mesh refinement (smallest cells).
"""

import vtk
import numpy as np
import math

def apply_magnification(renderer, active_actors, percentile_threshold=0.1, min_cells_required=10):
    """
    Apply view magnification focusing on regions with highest mesh refinement.
    
    Parameters:
    -----------
    renderer : vtkRenderer
        The VTK renderer to apply magnification to
    active_actors : list
        List of vtkActor objects representing active blocks
    percentile_threshold : float, optional
        Percentile threshold for cell size (default: 0.1 means focus on smallest 10% of cells)
    min_cells_required : int, optional
        Minimum number of cells required to consider a region (default: 10)
    """
    print("Analyzing mesh refinement patterns...")
    # Collect all cell sizes and their centers
    all_cell_sizes = []
    all_cell_centers = []
    
    for actor in active_actors:
        if not hasattr(actor, 'GetMapper') or not actor.GetMapper() or not actor.GetMapper().GetInput():
            continue
            
        data = actor.GetMapper().GetInput()
        
        # Skip if no cells
        if not hasattr(data, 'GetNumberOfCells') or data.GetNumberOfCells() == 0:
            continue
        
        # Process each cell to find its size (volume or area)
        num_cells = data.GetNumberOfCells()
        for i in range(num_cells):
            cell = data.GetCell(i)
            if not cell:
                continue
                
            # Calculate cell size based on type
            cell_size = 0
            cell_center = [0, 0, 0]
            
            # Use appropriate size metric based on cell type
            if cell.GetCellType() in [vtk.VTK_TETRA, vtk.VTK_VOXEL, vtk.VTK_HEXAHEDRON]:
                # 3D cells - use volume
                cell_size = cell.GetVolume()
            elif cell.GetCellType() in [vtk.VTK_TRIANGLE, vtk.VTK_QUAD]:
                # 2D cells - manually compute approximate area
                if cell.GetCellType() == vtk.VTK_TRIANGLE:
                    # For triangle, use triangle area formula with cross product
                    points = cell.GetPoints()
                    p0 = points.GetPoint(0)
                    p1 = points.GetPoint(1)
                    p2 = points.GetPoint(2)
                    
                    # Vectors for two sides
                    v1 = [p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]]
                    v2 = [p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]]
                    
                    # Cross product magnitude is twice the area
                    cross = [
                        v1[1]*v2[2] - v1[2]*v2[1],
                        v1[2]*v2[0] - v1[0]*v2[2],
                        v1[0]*v2[1] - v1[1]*v2[0]
                    ]
                    
                    cell_size = 0.5 * math.sqrt(cross[0]**2 + cross[1]**2 + cross[2]**2)
                else:  # QUAD
                    # For quad, use bounds to approximate area
                    bounds = [0] * 6
                    cell.GetBounds(bounds)
                    # Use the projected area on major planes
                    xy_area = (bounds[1] - bounds[0]) * (bounds[3] - bounds[2])
                    yz_area = (bounds[3] - bounds[2]) * (bounds[5] - bounds[4])
                    xz_area = (bounds[1] - bounds[0]) * (bounds[5] - bounds[4])
                    cell_size = max(xy_area, yz_area, xz_area)
            elif cell.GetCellType() in [vtk.VTK_LINE]:
                # 1D cells - use length
                p1 = cell.GetPoints().GetPoint(0)
                p2 = cell.GetPoints().GetPoint(1)
                dx = p2[0] - p1[0]
                dy = p2[1] - p1[1]
                dz = p2[2] - p1[2]
                cell_size = math.sqrt(dx*dx + dy*dy + dz*dz)
            else:
                # For other cell types, use a general approach
                bounds = [0] * 6
                cell.GetBounds(bounds)
                cell_size = ((bounds[1] - bounds[0]) * 
                            (bounds[3] - bounds[2]) * 
                            (bounds[5] - bounds[4]))
                
                if cell_size <= 0:
                    # Handle flat/degenerate cells by using max dimension
                    cell_size = max(bounds[1] - bounds[0],
                                   bounds[3] - bounds[2],
                                   bounds[5] - bounds[4])
            
            # Only include cells with valid size
            if cell_size > 0:
                # Get cell center
                points = cell.GetPoints()
                centroid = [0, 0, 0]
                for j in range(points.GetNumberOfPoints()):
                    point = points.GetPoint(j)
                    centroid[0] += point[0]
                    centroid[1] += point[1]
                    centroid[2] += point[2]
                    
                if points.GetNumberOfPoints() > 0:
                    centroid[0] /= points.GetNumberOfPoints()
                    centroid[1] /= points.GetNumberOfPoints()
                    centroid[2] /= points.GetNumberOfPoints()
                    
                    all_cell_sizes.append(cell_size)
                    all_cell_centers.append(centroid)
    
    if not all_cell_sizes or len(all_cell_sizes) < min_cells_required:
        print("Not enough valid cells found for mesh refinement magnification")
        return
    
    # Find cell size threshold for smallest cells (highest refinement)
    threshold = np.percentile(all_cell_sizes, percentile_threshold * 100)
    
    # Filter cells smaller than threshold (highest refinement)
    refined_centers = [
        all_cell_centers[i] 
        for i in range(len(all_cell_sizes)) 
        if all_cell_sizes[i] <= threshold
    ]
    
    if not refined_centers or len(refined_centers) < min_cells_required:
        print(f"Not enough cells below size threshold ({threshold})")
        return
    
    # Calculate bounds of refined region
    x_vals = [p[0] for p in refined_centers]
    y_vals = [p[1] for p in refined_centers]
    z_vals = [p[2] for p in refined_centers]
    
    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    z_min, z_max = min(z_vals), max(z_vals)
    
    # Calculate center of refined region
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
    
    # Calculate full range vs. refined range
    all_x = [c[0] for c in all_cell_centers]
    all_y = [c[1] for c in all_cell_centers]
    all_z = [c[2] for c in all_cell_centers]
    
    full_range = max(
        max(all_x) - min(all_x),
        max(all_y) - min(all_y),
        max(all_z) - min(all_z)
    )
    
    refined_range = max(
        x_max - x_min,
        y_max - y_min,
        z_max - z_min
    )
    
    # Calculate zoom factor - emphasize the refined region
    zoom_factor = 1.0
    if refined_range > 0:
        zoom_factor = full_range / refined_range * 0.8  # 80% to avoid cutting off edges
    #zoom_factor = min(max(zoom_factor, 1.0), 3.0)  # Limit to reasonable range
    
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
    
    print(f"Applied mesh refinement magnification with zoom factor: {zoom_factor:.2f}")
    print(f"Focus center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    print(f"Cell size threshold: {threshold:.6f}")
    print(f"Number of cells in refined region: {len(refined_centers)}")