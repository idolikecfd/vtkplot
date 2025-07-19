import argparse
import os.path
import sys
from . import magnification
from .core import process_tecplot_file


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
        original_apply = magnification.statistical.apply_magnification
        magnification.statistical.apply_magnification = lambda renderer, active_actors: original_apply(
            renderer, active_actors, percentile_cutoff=args.percentile_cutoff)
    
    elif args.magnification == 'variability':
        original_apply = magnification.variability.apply_magnification
        magnification.variability.apply_magnification = lambda renderer, active_actors: original_apply(
            renderer, active_actors, scalar_name=args.scalar_name, gradient_percentile=args.gradient_percentile)
    
    elif args.magnification == 'field_range':
        original_apply = magnification.field_range.apply_magnification
        magnification.field_range.apply_magnification = lambda renderer, active_actors: original_apply(
            renderer, active_actors, scalar_name=args.scalar_name, 
            min_percentile=args.min_percentile, max_percentile=args.max_percentile)
    
    elif args.magnification == 'bounding_box':
        original_apply = magnification.bounding_box.apply_magnification
        magnification.bounding_box.apply_magnification = lambda renderer, active_actors: original_apply(
            renderer, active_actors, padding_factor=args.padding_factor)
            
    elif args.magnification == 'mesh_refinement':
        # For mesh_refinement, check if we have the attributes (from subparser) or use defaults
        original_apply = magnification.mesh_refinement.apply_magnification
        if hasattr(args, 'percentile_threshold'):
            # User specified mesh_refinement subcommand with options
            magnification.mesh_refinement.apply_magnification = lambda renderer, active_actors: original_apply(
                renderer, active_actors, percentile_threshold=args.percentile_threshold, 
                min_cells_required=args.min_cells)
        else:
            # Default case - no subcommand specified, use default parameters
            magnification.mesh_refinement.apply_magnification = lambda renderer, active_actors: original_apply(
                renderer, active_actors, percentile_threshold=0.1, min_cells_required=10)
    
    # Process the Tecplot file with the specified magnification method
    process_tecplot_file(args.layout_file, args.magnification)


if __name__ == '__main__':
    main()