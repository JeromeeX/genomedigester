#!/usr/bin/env python3
# genome_analysis.py

import re

def read_fasta(file_path):
    """
    Reads a FASTA file and returns a dictionary where keys are the headers (chromosomes or enzyme names)
    and values are the corresponding sequences.
    """
    with open(file_path, 'r') as file:
        sequences = {}
        header = None # Initiating header to None
        sequence = []

        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences[header] = ''.join(sequence)
                header = line[1:] # Remove '>' character
                sequence = []
            else:
                sequence.append(line)
        
        if header: # Add the last sequence to the dictionary
            sequences[header] = ''.join(sequence)

    return sequences

def restriction_sites(file_path):
    """
    Reads the enzyme file and returns a dictionary where keys are the enzymes and values are a list
    consisting cutsite and the nick.
    """
    enzymes = read_fasta(file_path)
    for key in enzymes:
        nick = enzymes[key].find("|")
        enzymes[key] = [enzymes[key].replace("|", "").replace("N", "[ACTG]").replace("Y", "[CT]").replace("R", "[AG]"), nick] # Replace all the placeholders with the corresponding regex
    return enzymes

def find_restriction_sites(genome, enzyme):
    """
    Finds restriction sites in the genome. Returns a dictionary where keys are headers (chromosomes)
    and values are lists of positions of the restriction site.
    """
    positions = {}
    restriction_site = restriction_sites("test_enzyme.fa")[enzyme]

    for chr, seq in genome.items():
        positions[chr] = [0]+ [restriction_site[1] + m.start() for m in re.finditer(restriction_site[0], seq)] # Find all the restriction sites and add the nick to the position, add the first position of the chromosome
        if positions[chr][-1] != len(seq): # Add the last position of the chromosome
            positions[chr].append(len(seq))

    return positions

def calculate_distances(sites):
    """
    Calculates the distances between the restriction sites. Returns a dictionary where keys are
    headers (chromosomes) and values are lists of distances between the restriction sites.
    """
    distances = {}
    for chr, site in sites.items():
        distances[chr] = [site[i+1] - site[i] for i in range(len(site)-1)]
    return distances

def main():
    genome_file = "test_genome.fa"
    enzyme_file = "test_enzyme.fa"

    # Read genome and restriction enzyme sites
    genome = read_fasta(genome_file)
    enzymes = restriction_sites(enzyme_file)
    
    # Create a dictionary to store the restriction sites for the selected enzymes
    selected_enzymes = input("Enter the names of the enzymes separated by commas (e.g., EcoRI, HindIII): ").split(',')
    
    # Find restriction sites for the selected enzymes and generate the selected_restriction_sites dictionary
    from collections import defaultdict
    selected_restriction_sites = defaultdict(set) # Use defaultdict to avoid key errors and store unique restriction sites
    for enzyme_name in selected_enzymes:
        enzyme_name = enzyme_name.strip() # Remove whitespaces
        if enzyme_name in enzymes:
            enzyme_sites = find_restriction_sites(genome, enzyme_name)
            for chr, sites in enzyme_sites.items(): # Add the restriction sites to the selected_restriction_sites dictionary
                selected_restriction_sites[chr].extend(sites)
        else:
            print("{enzyme} not found!".format(enzyme=enzyme_name))
    
    for chr in selected_restriction_sites: # Sort the restriction sites and convert the set to sorted list
        selected_restriction_sites[chr] = sorted(list(selected_restriction_sites[chr]))
    
    # Calculate the distances between the restriction sites in each chromosome and merge together
    distances_in_chr = calculate_distances(selected_restriction_sites)
    distances = []
    for chr, distance in distances_in_chr.items():
        distances.extend(distance)
    
    return distances # Return the distances list to be used in the plot

if __name__ == "__main__":
    main()