import argparse
import time
import sys
import os
import csv
import re
from operator import itemgetter


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True,
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int,
                        default=50, help="Minimum gene length to consider")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int,
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes (shine box not included).")
    parser.add_argument('-p', dest='predicted_genes_file', type=str,
                        default=os.curdir + os.sep + "predict_genes.csv",
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=str,
                        default=os.curdir + os.sep + "genes.fna",
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file):
    """Extract the complete genome sequence as a single string
    """
    with open(fasta_file, 'r') as file:
        sequence = ''
        for line in file:
            if not str(line).startswith('>'):
                sequence += str(line).strip()
    sequence = sequence.upper()
    sequence = sequence.replace('U', 'T')
    return sequence


def find_start(start_regex, sequence, start, stop):
    """Find the start codon
    """
    match = start_regex.search(sequence, start, stop)
    if match is None:
        return None
    else:
        return match.start(0)


def find_stop(stop_regex, sequence, start):
    """Find the stop codon
    """
    match_reader = stop_regex.finditer(sequence, start)
    for match in match_reader:
        if match is None:
            return None
        else:
            # Est-ce que l'on est dans le cadre de lecture ?
            stop = match.start(0)
            if (stop-start) % 3 == 0:
                return stop


def has_shine_dalgarno(shine_regex, sequence, start,
                       max_shine_dalgarno_distance):
    """Find a shine dalgarno motif before the start codon
    """
    start_search_sd = start - max_shine_dalgarno_distance
    if start_search_sd < 0:
        start_search_sd = 0

    match = shine_regex.search(sequence, start_search_sd, start)
    if match is None:
        return False
    else:
        if (match.end(0) + 6) >= start:
            return False
        else:
            return True


def predict_genes(sequence, start_regex, stop_regex, shine_regex,
                  min_gene_len, max_shine_dalgarno_distance, min_gap):
    """Predict most probable genes
    """
    current_pos = 0
    genes_list = []

    while len(sequence) - current_pos >= min_gap:
        new_gene = False
        current_pos = find_start(start_regex, sequence, current_pos,
                                 len(sequence))

        if current_pos is not None:
            stop = find_stop(stop_regex, sequence, current_pos)
            if stop is not None:
                gene_len = stop - current_pos
                if gene_len >= min_gene_len:
                    if has_shine_dalgarno(shine_regex, sequence, current_pos,
                                          max_shine_dalgarno_distance):
                        new_gene = True  # On a trouvé un gène !

        if new_gene:
            genes_list.append([current_pos + 1, stop + 3])
            current_pos = stop + 2 + min_gap
        else:
            current_pos += 1

    return genes_list


def write_genes_pos(predicted_genes_file, probable_genes):
    """Write list of gene positions
    """
    try:
        with open(predicted_genes_file, "wt") as genes_file:
            predict_genes_writer = csv.writer(genes_file, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_genes(fasta_file, sequence, probable_genes, sequence_rc,
                probable_genes_comp):
    """Write gene sequence in fasta format
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i, gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i+1, os.linesep,
                    fill(sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i+1
            for j, gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(kmer):
    """Get the reverse complement"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in kmer[::-1]])


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """

    parse = False

    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    #start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    #stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    # Shine AGGAGGUAA
    # AGGA ou GGAGG
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')

    # Arguments
    args = get_arguments()
    genome_file = args.genome_file
    min_gene_len = args.min_gene_len
    max_shine_dalgarno_distance = args.max_shine_dalgarno_distance
    min_gap = args.min_gap
    predicted_genes_file = args.predicted_genes_file
    fasta_file = args.fasta_file

    # Sens 5' -> 3'
    sequence = read_fasta(genome_file)
    genes_53 = predict_genes(sequence, start_regex, stop_regex, shine_regex,
                             min_gene_len,
                             max_shine_dalgarno_distance,
                             min_gap)
    print(f"Nb de gènes 5' -> 3' : {len(genes_53)}")

    # Sens 3' -> 5'
    sequence_rc = reverse_complement(sequence)
    genes_35 = predict_genes(sequence_rc, start_regex, stop_regex, shine_regex,
                             min_gene_len,
                             max_shine_dalgarno_distance,
                             min_gap)
    print(f"Nb de gènes 3' -> 5' : {len(genes_35)}")
    genes_35 = [[len(sequence) - e + 1, len(sequence) - s + 1]
                for s, e in genes_35][::-1]

    # Ecriture des résultats
    write_genes_pos(predicted_genes_file,
                    sorted(genes_53 + genes_35, key=itemgetter(0)))
    write_genes(fasta_file, sequence, genes_53, sequence_rc, genes_35)


if __name__ == '__main__':
    start_time = time.time()
    main()
    end_time = time.time() - start_time
    print(f'Genes prediction done in {end_time // 60:.0f} min'
          f' {end_time % 60:.2f} s.')
