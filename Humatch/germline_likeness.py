import os
import requests
import numpy as np
from Humatch.utils import (
    get_ordered_AA_one_letter_codes,
    get_CDR_loop_indices,
    get_indices_of_selected_imgt_positions_in_canonical_numbering
    )

# gl arrays added to compiled env package_data
HUMATCH_CODE_DIR = os.path.dirname(os.path.abspath(__file__))
# set HUMATCH_WEIGHTS_DIR to keep trained_models/ and germline_likeness_lookup_arrays/ elsewhere
HUMATCH_WEIGHTS_DIR = os.environ.get("HUMATCH_WEIGHTS_DIR") or HUMATCH_CODE_DIR
GL_DIR = os.path.join(HUMATCH_WEIGHTS_DIR, "germline_likeness_lookup_arrays")
if not os.path.exists(GL_DIR):
    os.makedirs(GL_DIR)
vgenes = [f"hv{i}" for i in range(1, 8)] + [f"lv{i}" for i in range(1, 11)] + [f"kv{i}" for i in range(1, 8)]
gl_urls = {gene: f"https://zenodo.org/records/13764771/files/{gene}.npy?download=1" for gene in vgenes}


def mutate_seq_to_match_germline_likeness(seq, target_gene, target_score, allow_CDR_mutations=False,
                                          fixed_imgt_positions=[], germline_likeness_lookup_arrays_dir=GL_DIR):
    '''
    Iteratively mutate a sequence to achieve a target germline likeness

    :param seq: str, sequence (padded with "-" for missing positions)
    :param germline_likeness_lookup_arrays_dir: str, path to directory containing the previously
        calculated position AA frequencies
    :param target_gene: str, target gene e.g. "hv1", "lv3", "kv6" etc.
    :param target_score: float, target germline likeness score
    :param allow_CDR_mutations: bool, if CDR mutations are allowed
    :param fixed_imgt_positions: list, list of IMGT positions to exclude from mutation
    :return: str, mutated sequence
    '''
    # get starting germline likeness score
    germline_score = get_normalised_germline_likeness_score(seq, target_gene, germline_likeness_lookup_arrays_dir)
    while germline_score < target_score:
        # make top single mutation to most improve germline likeness
        new_seq = make_top_N_most_observed_germline_mutations(seq, target_gene, 1, germline_likeness_lookup_arrays_dir=germline_likeness_lookup_arrays_dir,
                                                              allow_CDR_mutations=allow_CDR_mutations, fixed_imgt_positions=fixed_imgt_positions)
        # break if no more differing positions found between sequence and germline
        if new_seq == seq:
            break
        # update germline score
        seq = new_seq
        germline_score = get_normalised_germline_likeness_score(seq, target_gene, germline_likeness_lookup_arrays_dir)
    return seq


def load_observed_position_AA_freqs(target_gene, germline_likeness_lookup_arrays_dir=GL_DIR):
    '''
    Load the observed position AA frequencies for a target gene

    :param germline_likeness_lookup_arrays_dir: str, path to directory
        containing the previously calculated position AA frequencies
    :param target_gene: str, target gene e.g. "hv1", "lv3", "kv6" etc.
    :return: np.array, observed position AA frequencies, shape (200, 20)
    '''
    lookup_arr = os.path.join(germline_likeness_lookup_arrays_dir, f"{target_gene}.npy")
    if not os.path.exists(lookup_arr):
        print(f"Downloading germline lookup arrays for all V-genes from https://zenodo.org/records/13764771... ", end="")
        for gene, url in gl_urls.items():
            r = requests.get(url)
            with open(os.path.join(germline_likeness_lookup_arrays_dir, f"{gene}.npy"), 'wb') as f:
                f.write(r.content)
        print("done")
    return np.load(lookup_arr)


def get_list_of_occurence_freqs_for_seq_based_on_gene_arr(seq, arr):
    '''
    Get the observed position AA frequencies for a given sequence

    :param seq: str, sequence (padded with "-" for missing positions)
    :param arr: np.array, observed position AA frequencies, shape (200, 20)
    :return: np.array, observed position AA frequencies for the given sequence, shape (200,)
    '''
    # covert str seq to list of index positions (A,C,D,... -> 0,1,2,...)
    AA_indexes = [get_ordered_AA_one_letter_codes().index(i) for i in seq]
    # add pad dim to arr filled with zeros - to account for "-" in seqs, new shape (200, 21)
    arr = np.pad(arr, ((0,0), (0,1)), 'constant', constant_values=(0))
    # get values from arr at AA_indexes
    freqs = arr[np.arange(len(arr)), AA_indexes]
    return freqs


def get_most_common_germline_seq(arr):
    '''
    Get the most `common' germline sequence based on the observed position AA frequencies
        i.e. take the most common AA at each position
        When no AA is observed at a position (e.g. large insertions in CDRH3), idx 0 is taken

    :param arr: np.array, observed position AA frequencies, shape (200, 20)
    :return: str, most common germline sequence
    '''
    # get indexes of max values for arr e.g. 200x20 --> 200x1
    indices_of_most_common_AAs = np.argmax(arr, axis=1)
    # get AA one letter codes for max_idx
    AA_codes = [get_ordered_AA_one_letter_codes()[i] for i in indices_of_most_common_AAs]
    return "".join(AA_codes)


def get_normalised_germline_likeness_score(seq, target_gene, germline_likeness_lookup_arrays_dir=GL_DIR):
    '''
    We define the germline likeness of a sequence as the average observed position AA frequency
        We divide by the padded length of the sequence to normalise the germline likeness as long/short
        CDRH3s can have outsized effects when it should not as it is (largely) not V-gene encoded

    :param seq: str, sequence (padded with "-" for missing positions)
    :param germline_likeness_lookup_arrays_dir: str, path to directory containing the previously
        calculated position AA frequencies
    :param target_gene: str, target gene e.g. "hv1", "lv3", "kv6" etc.
    :return: float, normalised germline likeness of the sequence
    '''
    arr = load_observed_position_AA_freqs(target_gene, germline_likeness_lookup_arrays_dir)
    freqs = get_list_of_occurence_freqs_for_seq_based_on_gene_arr(seq, arr)
    return np.sum(freqs) / len(seq)


def get_indices_where_two_strs_do_not_match(str1, str2, pad_chars=["-", "X", "*"]):
    return [i for i in range(len(str1)) if (str1[i] != str2[i]) and (str1[i] not in pad_chars) and (str2[i] not in pad_chars)]


def get_ranked_indices_to_mutate(seq, arr, allow_CDR_mutations=False, fixed_imgt_positions=[]):
    '''
    Get the indices that look least like the germline sequence

    :param seq: str, sequence (padded with "-" for missing positions)
    :param arr: np.array, observed position AA frequencies, shape (200, 20)
    :param allow_CDR_mutations: bool, if CDR mutations are allowed
    :param fixed_imgt_positions: list, list of IMGT positions to exclude from mutation
    :return: list, ranked indices to mutate
    '''
    # rank indices based on difference in occurence freqs - interested in the most differing positions
    germline_seq = get_most_common_germline_seq(arr)
    indices_where_seq_and_germline_diff = get_indices_where_two_strs_do_not_match(seq, germline_seq)
    occurence_freqs_of_starting_seq = get_list_of_occurence_freqs_for_seq_based_on_gene_arr(seq, arr)
    occurence_freqs_of_germline_seq = get_list_of_occurence_freqs_for_seq_based_on_gene_arr(germline_seq, arr)
    difference_in_freqs = occurence_freqs_of_germline_seq - occurence_freqs_of_starting_seq
    ranked_indices = sorted(indices_where_seq_and_germline_diff, key=lambda x: difference_in_freqs[x], reverse=True)
    # remove CDR loop indices and others from consideration if not allowed
    fixed_indices = []
    if not allow_CDR_mutations:
        fixed_indices += get_CDR_loop_indices()
    if fixed_imgt_positions:
        fixed_indices += get_indices_of_selected_imgt_positions_in_canonical_numbering(fixed_imgt_positions)
    ranked_indices = [i for i in ranked_indices if i not in fixed_indices]
    return ranked_indices


def make_top_N_most_observed_germline_mutations(seq, target_gene, N, allow_CDR_mutations=False,
                                                fixed_imgt_positions=[], germline_likeness_lookup_arrays_dir=GL_DIR):
    '''
    Make the top N most observed germline mutations to a sequence (we use it for N=1)

    :param seq: str, sequence (padded with "-" for missing positions)
    :param germline_likeness_lookup_arrays_dir: str, path to directory containing the previously
        calculated position AA frequencies
    :param target_gene: str, target gene e.g. "hv1", "lv3", "kv6" etc.
    :param N: int, number of mutations to make
    :param allow_CDR_mutations: bool, if CDR mutations are allowed
    :param fixed_imgt_positions: list, list of IMGT positions to exclude from mutation
    :return: str, mutated sequence
    '''
    arr = load_observed_position_AA_freqs(target_gene, germline_likeness_lookup_arrays_dir)
    germline_seq = get_most_common_germline_seq(arr)
    ranked_indices = get_ranked_indices_to_mutate(seq, arr, allow_CDR_mutations=allow_CDR_mutations,
                                                  fixed_imgt_positions=fixed_imgt_positions)
    # avoid infinite loop if no differing positions found between sequence and germline
    if len(ranked_indices) == 0:
        print("Warning: No differing positions found between sequence and germline")
        print("No mutations made")
        return seq
    elif N > len(ranked_indices):
        print(f"Warning: N ({N}) is greater than the number of differing positions ({len(ranked_indices)})")
        print(f"Mutating maximum number of differing positions ({len(ranked_indices)}) only")
        N = len(ranked_indices)
    else:
        pass
    # make top N mutations
    for i in range(N):
        seq = seq[:ranked_indices[i]] + germline_seq[ranked_indices[i]] + seq[ranked_indices[i]+1:]
    return seq
