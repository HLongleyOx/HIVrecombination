#!/usr/bin/env python3

"""
SLIM Output to FASTA Converter

This script converts SLIM simulation output files (.txt) to FASTA format,
handling both segregating and fixed mutations.

Usage:
    python slim_to_fasta.py <generation> <replicate> <rho_value> [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR] [--ref-fasta REF_FASTA]

Arguments:
    generation:  Generation number from SLIM simulation
    replicate:   Replicate number of the simulation
    rho_value:   Rho value used in the simulation

Options:
    --input-dir:  Directory containing SLIM output files (default: ./slim_output)
    --output-dir: Directory for output FASTA files (default: ./fasta_output)
    --ref-fasta:  Path to reference FASTA file (default: ./reference/HXB2.fasta)
    
Example:
    python slim_to_fasta.py 1000 1 0.5 --input-dir ./my_slim_data --output-dir ./my_fasta --ref-fasta ./refs/HXB2.fasta

Directory Structure:
    $SLIM_BASE_DIR/
    ├── reference/
    │   └── HXB2.fasta
    ├── simulations/
    │   └── simulation_rep{replicate}_rho{rho}/
    │       ├── gen_{generation}_seg.txt
    │       └── gen_{generation}_fixed.txt
    └── sequences/
        └── output_{replicate}_rho{rho}_{generation}.fasta
"""

import os
import sys
from typing import List, Tuple, Optional

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def parse_segregating_file(filepath: str) -> Tuple[List[Tuple[int, str, int]], List[Tuple[int, List[str]]]]:
    """
    Parse SLIM output file containing segregating mutations.

    Args:
        filepath: Path to the SLIM output file

    Returns:
        Tuple containing:
        - List of (mutation_id, nucleotide, position) tuples
        - List of (genome_id, mutations) tuples

    Raises:
        FileNotFoundError: If the input file doesn't exist
        ValueError: If file format is incorrect
    """
    try:
        with open(filepath, 'r') as file:
            lines = file.readlines()

        # Find section markers
        try:
            mutations_start = lines.index('Mutations:\n') + 1
            genomes_start = lines.index('Genomes:\n')
        except ValueError:
            raise ValueError("File format error: Missing required section headers")

        # Parse mutations
        mutations_data = []
        for line in lines[mutations_start:genomes_start]:
            fields = line.strip().split()
            mutations_data.append((
                int(fields[0]),      # mutation_id
                fields[9],           # nucleotide
                int(fields[3])       # position
            ))

        # Parse genomes
        genomes_data = []
        for line in lines[genomes_start + 1:]:
            fields = line.strip().split()
            genome_id = int(fields[0].split(':')[1])
            genomes_data.append((genome_id, fields[2:]))

        return mutations_data, genomes_data

    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find file: {filepath}")


def parse_fixed_file(filepath: str) -> List[Tuple[str, int]]:
    """
    Parse SLIM output file containing fixed mutations.

    Args:
        filepath: Path to the SLIM output file

    Returns:
        List of (nucleotide, position) tuples for fixed mutations

    Raises:
        FileNotFoundError: If the input file doesn't exist
        ValueError: If file format is incorrect
    """
    try:
        with open(filepath, 'r') as file:
            lines = file.readlines()

        try:
            mutations_start = lines.index('Mutations:\n') + 1
        except ValueError:
            raise ValueError("File format error: Missing 'Mutations:' header")

        return [(line.strip().split()[9], int(line.strip().split()[3]))
                for line in lines[mutations_start:]]

    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find file: {filepath}")


def read_fasta(filepath: str) -> List[str]:
    """
    Read sequences from a FASTA file.

    Args:
        filepath: Path to the FASTA file

    Returns:
        List of sequences

    Raises:
        FileNotFoundError: If the input file doesn't exist
    """
    try:
        with open(filepath) as handle:
            return [record for _, record in SeqIO.FastaIO.SimpleFastaParser(handle)]
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find FASTA file: {filepath}")


def write_sequences_to_fasta(sequences: List[List[str]], 
                           output_file: str, 
                           time: int) -> None:
    """
    Write sequences to a FASTA file.

    Args:
        sequences: List of sequences to write
        output_file: Path to output FASTA file
        time: Time point for sequence identification

    Raises:
        IOError: If unable to write to output file
    """
    try:
        with open(output_file, "w") as f:
            for idx, seq in enumerate(sequences):
                seq_obj = Seq(seq[0])
                seq_record = SeqRecord(
                    seq_obj,
                    id=f"sequence_{time}_{idx}",
                    description=""
                )
                SeqIO.write(seq_record, f, "fasta")
    except IOError as e:
        raise IOError(f"Error writing to {output_file}: {str(e)}")


def process_sequence(sequence: str,
                    mutations_data: List[Tuple[int, str, int]],
                    genomes_data: List[Tuple[int, List[str]]],
                    fixed_data: List[Tuple[str, int]]) -> List[List[str]]:
    """
    Process the ancestral sequence with mutations to generate new sequences.

    Args:
        sequence: Original ancestral sequence
        mutations_data: List of segregating mutations
        genomes_data: List of genome data
        fixed_data: List of fixed mutations

    Returns:
        List of processed sequences
    """
    # Create ancestral sequence list
    ancestral_sequence = [(base, i) for i, base in enumerate(sequence)]
    
    genomes_all = []
    
    for genome_data in genomes_data:
        segregating_muts = genome_data[1]
        genome = ""

        # Convert mutation IDs to integers
        muts_numeric = [int(mut) for mut in segregating_muts]
        
        # Get segregating alleles
        segregating_allele = [
            [nuc, pos] for mut_id, nuc, pos in mutations_data 
            if mut_id in muts_numeric
        ]

        # Process each position
        for i, (anc_base, _) in enumerate(ancestral_sequence):
            # Check segregating mutations
            seg_mutations = [nuc for nuc, pos in segregating_allele if i == pos]
            # Check fixed mutations
            fix_mutations = [nuc for nuc, pos in fixed_data if i == pos]
            
            if seg_mutations:
                genome += seg_mutations[0]
            elif fix_mutations:
                genome += fix_mutations[0]
            else:
                genome += anc_base

        genomes_all.append([genome])

    return genomes_all


def parse_args():
    """Parse command line arguments."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Convert SLIM output to FASTA format")
    parser.add_argument("generation", type=int, help="Generation number from SLIM simulation")
    parser.add_argument("replicate", type=int, help="Replicate number of the simulation")
    parser.add_argument("rho_value", type=float, help="Rho value used in the simulation")
    parser.add_argument("--input-dir", default="./slim_output", 
                       help="Directory containing SLIM output files")
    parser.add_argument("--output-dir", default="./fasta_output",
                       help="Directory for output FASTA files")
    parser.add_argument("--ref-fasta", default="./reference/HXB2.fasta",
                       help="Path to reference FASTA file")
    
    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_args()

    # Convert arguments to proper types and names
    generation = str(args.generation)
    replicate = str(args.replicate)
    rho = str(args.rho_value)

    # Set up and validate paths
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    ref_fasta = os.path.abspath(args.ref_fasta)

    # Ensure required directories exist
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(ref_fasta), exist_ok=True)
    
    try:
        # Read input files
        sequence = read_fasta(ref_fasta)[0]
        
        # Define input and output files
        seg_file = os.path.join(input_dir, f'gen_{generation}_seg.txt')
        fixed_file = os.path.join(input_dir, f'gen_{generation}_fixed.txt')
        
        # Parse input files
        mutations_data, genomes_data = parse_segregating_file(seg_file)
        fixed_data = parse_fixed_file(fixed_file)

        # Process sequences
        genomes_all = process_sequence(sequence, mutations_data, genomes_data, fixed_data)

        # Calculate time and set up output path
        time = int(generation) - 50000
        output_file = os.path.join(
            output_dir,
            f"output_{replicate}_rho{rho}_{generation}.fasta"
        )
        
        write_sequences_to_fasta(genomes_all, output_file, time)
        print(f"Successfully converted sequences to: {output_file}")

    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()