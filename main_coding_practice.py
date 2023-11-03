import pandas as pd
import ahocorasick as acs
from bokeh.plotting import figure, show, output_file
import gzip
import random

# Read the peptide column from psm.tsv.gz
with gzip.open(r"psm.tsv.gz", "rt") as fin:
    peptides_df = pd.read_csv(fin, sep='\t', usecols=['Peptide'])

# Aho-Corasick automaton for peptides
peptides = acs.Automaton()
for peptide in peptides_df['Peptide']:
    peptides.add_word(peptide, peptide)
peptides.make_automaton()

# dictionary to store peptide-protein mappings
matched_peptides = {}

# Read in proteins from UniProt_Human.fasta.gz
with gzip.open(r"UniProt_Human.fasta.gz", "rt") as fin:
    protein_sequences = fin.read()

    # Iterate over matches
    for end_index, peptide in peptides.iter(protein_sequences):
        start_index = end_index - len(peptide) + 1

        # Extract protein entry from protein sequences
        protein_entry = protein_sequences.rfind('>', 0, start_index)
        protein_entry = protein_sequences[protein_entry:end_index].strip()

        # Extract protein ID and sequence
        lines = protein_entry.split('\n')
        protein_id = lines[0].split()[0]
        protein_sequence = ''.join(lines[1:])

        if peptide in protein_sequence:
            if peptide not in matched_peptides:
                matched_peptides[peptide] = []
            matched_peptides[peptide].append(protein_id)

# Print peptide-protein mappings
for peptide, proteins in matched_peptides.items():
    print(f"Peptide: {peptide}, Proteins: {', '.join(proteins)}")

import pandas as pd
import ahocorasick as acs
import gzip

# Read the peptide column from psm.tsv.gz
with gzip.open(r"psm.tsv.gz", "rt") as fin:
    peptides_df = pd.read_csv(fin, sep='\t', usecols=['Peptide'])

# Create Aho-Corasick automaton for peptides
peptides = acs.Automaton()
for peptide in peptides_df['Peptide']:
    peptides.add_word(peptide, peptide)
peptides.make_automaton()

# Create dictionary to store peptide-protein mappings
matched_peptides = {}

# Read in proteins from UniProt_Human.fasta.gz
with gzip.open(r"UniProt_Human.fasta.gz", "rt") as fin:
    protein_sequences = fin.read()

    # Iterate over the matches found by the automaton
    for end_index, peptide in peptides.iter(protein_sequences):
        start_index = end_index - len(peptide) + 1

        # Extract protein entry from the protein sequences
        protein_entry = protein_sequences.rfind('>', 0, start_index)
        protein_entry = protein_sequences[protein_entry:end_index].strip()

        # Extract protein ID and sequence
        lines = protein_entry.split('\n')
        protein_id = lines[0].split()[0]
        protein_sequence = ''.join(lines[1:])

        if peptide in protein_sequence:
            if peptide not in matched_peptides:
                matched_peptides[peptide] = []
            matched_peptides[peptide].append(protein_id)

# Print peptide-protein mappings
for peptide, proteins in matched_peptides.items():
   print(f"Peptide: {peptide}, Proteins: {', '.join(proteins)}")

# frequency of each identified protein
protein_counts = {}
for proteins in matched_peptides.values():
    for protein in proteins:
        if protein not in protein_counts:
            protein_counts[protein] = 0
        protein_counts[protein] += 1

# dataframe for the protein counts
protein_df = pd.DataFrame(list(protein_counts.items()), columns=['Protein', 'Frequency'])

# Sort the dataframe by frequency in descending order
protein_df = protein_df.sort_values(by='Frequency', ascending=False)

# Bokeh figure
p = figure(x_range=protein_df['Protein'], height=600, width=1000, title="Protein Frequency Coverage Map",
           toolbar_location=None, tools="")

# appearance of the bar graph
p.vbar(x='Protein', top='Frequency', width=0.9, source=protein_df, line_color='white')

# rotate x-axis labels
p.xaxis.major_label_orientation = 1.2

# Set axis labels and plot title
p.xaxis.axis_label = "Protein ID"
p.yaxis.axis_label = "Frequency"
p.title.align = "center"

# Show the plot
show(p)

# Create a dictionary to store protein UniProt ID as key and identified PSMs as values
output_dict = {}

# Iterate over matched peptides and proteins
for peptide, proteins in matched_peptides.items():
    for protein in proteins:
        # Check if the protein ID is already in the dictionary, if not, add it with an empty list
        if protein not in output_dict:
            output_dict[protein] = []
        # Add the current peptide to the list of identified PSMs for the protein
        output_dict[protein].append(peptide)

# Print the output dictionary
for protein_id, peptides in output_dict.items():
    print(f"protein uniprot_id: {protein_id}, identified psm(s): {', '.join(peptides)}")

# Function to reverse a sequence
def reverse_sequence(sequence):
    return sequence[::-1]

# Function to reverse a segment based on cut site (KR)
def reverse_segment(sequence, cut_site="KR"):
    segments = sequence.split(cut_site)
    reversed_segments = [segment[::-1] for segment in segments]
    return cut_site.join(reversed_segments)

# Function to randomly rearrange a sequence
def random_rearrange(sequence):
    sequence_list = list(sequence)
    random.shuffle(sequence_list)
    return ''.join(sequence_list)

# fasta file
input_file = "input.fasta"

# options: "A" for whole sequence reverse, "B" for segmented reverse, "C" for random rearrangement
option = "B"

output_sequences = []
with gzip.open(r"UniProt_Human.fasta.gz", "rt") as fin:
    current_sequence_id = ""
    current_sequence = ""
    for line in fin:
        if line.startswith(">"):
            if current_sequence_id:
                if option == "A":
                    reverse_seq = reverse_sequence(current_sequence)
                elif option == "B":
                    reverse_seq = reverse_segment(current_sequence, cut_site="KR")
                elif option == "C":
                    reverse_seq = random_rearrange(current_sequence)
                else:
                    print("invalid")
                    break

                # original sequence
                output_sequences.append(f">{current_sequence_id}\n{current_sequence}\n")
                # reverse sequence
                output_sequences.append(f">rev_{current_sequence_id}\n{reverse_seq}\n")

            current_sequence_id = line.strip()[1:]
            current_sequence = ""
        else:
            current_sequence += line.strip()

    # last sequence
    if current_sequence_id:
        if option == "A":
            reverse_seq = reverse_sequence(current_sequence)
        elif option == "B":
            reverse_seq = reverse_segment(current_sequence, cut_site="KR")
        elif option == "C":
            reverse_seq = random_rearrange(current_sequence)
        else:
            print("invalid")

        # original sequence
        output_sequences.append(f">{current_sequence_id}\n{current_sequence}\n")
        # reverse sequence
        output_sequences.append(f">rev_{current_sequence_id}\n{reverse_seq}\n")

# Print original sequence
print(f">{current_sequence_id}")
print(current_sequence)
# Print reverse sequence
print(f">rev_{current_sequence_id}")
print(reverse_seq)