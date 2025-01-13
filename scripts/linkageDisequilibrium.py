#!/usr/bin/env python3

"""
Linkage Disequilibrium Analysis Tool

This script calculates linkage disequilibrium for recombination rate analysis
by identifying diverse sites in genetic sequences.

Usage:
    python ld_analysis.py <sim_number> <rho_value> [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR]

Arguments:
    sim_number: Simulation number
    rho_value: Rho value for the analysis

Options:
    --input-dir:  Directory containing input FASTA files (default: ./sequences)
    --output-dir: Directory for output CSV files (default: ./linkage)

Example:
    python ld_analysis.py 1 0.5 --input-dir ./my_sequences --output-dir ./my_results

Author: Harriet Longley
"""

import os
import sys
import argparse
from collections import Counter, defaultdict
from typing import List, Tuple, Dict

import pandas as pd
from Bio import SeqIO


[... previous functions remain unchanged ...]


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate linkage disequilibrium from sequence data"
    )
    
    parser.add_argument(
        "sim_number",
        help="Simulation number"
    )
    parser.add_argument(
        "rho_value",
        help="Rho value for the analysis"
    )
    parser.add_argument(
        "--input-dir",
        default="./sequences",
        help="Directory containing input FASTA files"
    )
    parser.add_argument(
        "--output-dir",
        default="./linkage",
        help="Directory for output CSV files"
    )
    
    return parser.parse_args()


def main():
    """Main execution function."""
    # Parse command line arguments
    args = parse_args()
    
    # Set up paths
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Construct input and output file paths
    fasta_file = os.path.join(
        input_dir,
        f"sequences_rho{args.rho_value}_{args.sim_number}.fasta"
    )
    
    # Read sequences
    try:
        with open(fasta_file) as handle:
            sample_sequences = [
                [index, record] 
                for index, record in SeqIO.FastaIO.SimpleFastaParser(handle)
            ]
    except FileNotFoundError:
        print(f"Error: Could not find file {fasta_file}")
        sys.exit(1)
    
    # Analysis parameters
    dates = [f"_{x}_" for x in range(0, 601, 50)]
    
    # Find diverse sites
    sites_diverse, diversity = find_sites(sample_sequences, dates)
    
    # Process diverse sites
    sites_diverse['date_str'] = sites_diverse['date']
    sites_diverse = sites_diverse.sort_values(by=["site", "date"])
    
    # Get first occurrence of each site
    mask = sites_diverse['date'] == sites_diverse.groupby('site')['date'].transform('min')
    sites_diverse_first = (
        sites_diverse[mask]
        .rename(columns={'nucleotide': 'first.nucleotide', 'date': 'first.date'})
        .drop("date_str", axis=1)
    )
    sites_diverse = pd.merge(sites_diverse, sites_diverse_first, on="site")
    
    # Calculate linkage
    df_linkage = pd.DataFrame()
    
    for i, df in sites_diverse.iterrows():
        base_a = df["first.nucleotide"]
        site_a = df["site"]
        current_date = df["date"]
        
        # Get other sites at same timepoint
        df_othersites = (
            sites_diverse
            .iloc[(i + 1):]
            [sites_diverse["date"] == current_date]
        )
        
        for _, df_b in df_othersites.iterrows():
            base_b = df_b["first.nucleotide"]
            site_b = df_b["site"]
            
            linkage_output = calculate_linkage(
                sample_sequences, 
                current_date, 
                site_a, 
                site_b, 
                base_a,
                base_b
            )
            
            if linkage_output:
                df_linkage = pd.concat([df_linkage, linkage_output])
    
    # Add metadata
    df_linkage["ID"] = args.sim_number
    df_linkage["rho"] = args.rho_value
    
    # Save results
    output_file = os.path.join(
        output_dir,
        f"linkage_{args.sim_number}_rho{args.rho_value}.csv"
    )
    df_linkage.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    main()