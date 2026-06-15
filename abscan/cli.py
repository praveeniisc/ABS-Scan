import argparse
import logging
import sys
from .scanner import AlanineScanner
from .plotting import plot_results

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )

def main():
    parser = argparse.ArgumentParser(
        description="ABS-Scan: Alanine Binding Site Mutagenesis Scanning for Protein-Ligand Complexes."
    )
    
    parser.add_argument(
        '-f', '--pdb',
        action='store',
        required=True,
        dest='pdbfile',
        help='PDB file of the protein-ligand complex'
    )
    
    parser.add_argument(
        '-n', '--resno',
        action='store',
        required=True,
        dest='resno',
        help='Residue number of the HETATM ligand'
    )
    
    parser.add_argument(
        '-d', '--dist',
        action='store',
        default=4.5,
        type=float,
        dest='dist',
        help='Distance cutoff (Å) to define the pocket residues (default: 4.5)'
    )
    
    parser.add_argument(
        '-o', '--outdir',
        action='store',
        required=True,
        dest='outdir',
        help='Output directory to store structures and results'
    )
    
    parser.add_argument(
        '-s', '--sf',
        action='store',
        default='vina',
        dest='sf',
        choices=['vina', 'ad4', 'vinardo'],
        help='AutoDock Vina scoring function to use (choices: vina, ad4, vinardo; default: vina)'
    )
    
    args = parser.parse_args()
    
    setup_logging()
    logger = logging.getLogger("abscan.cli")
    
    logger.info("Starting ABS-Scan pipeline...")
    logger.info(f"Input PDB: {args.pdbfile}")
    logger.info(f"Ligand residue number: {args.resno}")
    logger.info(f"Pocket cutoff distance: {args.dist} Å")
    logger.info(f"Scoring function: {args.sf}")
    logger.info(f"Output directory: {args.outdir}")
    
    try:
        scanner = AlanineScanner(
            pdb_file=args.pdbfile,
            resno=args.resno,
            distance=args.dist,
            outdir=args.outdir,
            scoring_function=args.sf
        )
        csv_path = scanner.run()
        
        logger.info("Generating visualization plots...")
        plot_results(csv_path, args.outdir)
        
        logger.info(f"ABS-Scan finished successfully! Results are saved in {args.outdir}")
        
    except Exception as e:
        logger.error(f"ABS-Scan pipeline failed: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
