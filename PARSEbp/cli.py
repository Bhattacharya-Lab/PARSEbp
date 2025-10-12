import argparse
from .core import parsebp

def main():
    
    parser = argparse.ArgumentParser(
        prog="PARSEbp",
        description=(
            "PARSEbp (Pairwise Agreement-based RNA Scoring with Emphasis on Base Pairings)"
        ),
        epilog="Example: PARSEbp --pdb_dir Inputs --output score.txt --mode 1 --num_threads 50 --target_seq """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--mode", type=int, default=1,
        help="Scoring mode: 1 (default) uses base pair guidance extracted from decoys; "
             "0 disables base pair guidance."
    )
    parser.add_argument(
        "--num_threads", type=int, default=50,
        help="Number of threads to use for parallel computation."
    )
    parser.add_argument(
        "--pdb_dir", type=str,  default="./Inputs",
        help="Path to the directory containing the RNA 3D structural ensemble (PDB files)."
    )
    parser.add_argument(
        "--output", type=str, default="./score.txt",
        help="Path to the output file where scores will be saved."
    )
    parser.add_argument(
        "--target_seq", type=str,  default="",
        help="Target sequence to match and filter the pdb files in the specified directory.(default = "")"
    )

    args = parser.parse_args()

    p = parsebp()
    p.set_mode(args.mode)
    p.set_parallel_threads(args.num_threads)
    p.set_target_sequnece(args.target_seq)

    p.load_pdbs(args.pdb_dir)
    
    score = p.score()
    score.save(args.output)
    
