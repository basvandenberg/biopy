import itertools
import numpy

# amino acids (www.ebi.ac.uk/2can/biology/molecules_small_aatable.html)
# order is important! do not change
aa_unambiguous_alph = 'ARNDCEQGHILKMFPSTWYV'

aa_unambiguous_short = [
    'ala', 'arg', 'asn', 'asp', 'cys', 'glu', 'gln', 'gly', 'his', 'ile',
    'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val']

aa_unambiguous_name = [
    'alanine',       'arginine',  'asparagine', 'aspartic acid', 'cysteine',
    'glutamic acid', 'glutamine', 'glycine',    'histidine',     'isoleucine',
    'leucine',       'lysine',    'methionine', 'phenylalanine', 'proline',
    'serine',        'threonine', 'tryptophan', 'tyrosine',      'valine']

aa_ambiguous_alph = 'BJZX'
aa_ambiguous_short = ['asx', 'xle', 'xaa', 'glx']
aa_ambiguous_name = [
    'aspartic acid or asparagine', 'leucine or isoleucine',
    'unknown amino acid', 'glutamic acid or glutamine'
]

aa_special_alph = 'UO'
aa_special_short = ['sec', 'pyl']
aa_special_name = ['selenocysteine', 'pyrralysine']

aa_ter_alph = '*'
aa_ter_short = ['ter']
aa_ter_name = ['terminal']

aa_alph = aa_unambiguous_alph + aa_ambiguous_alph + aa_special_alph +\
        aa_ter_alph
aa_short = list(itertools.chain.from_iterable([aa_unambiguous_short,
        aa_ambiguous_short, aa_special_short, aa_ter_short]))
aa_name = list(itertools.chain.from_iterable([aa_unambiguous_name,
        aa_ambiguous_name, aa_special_name, aa_ter_name]))

# human mutations: key is from, values are occuring to mutations
non_zero_mutation_counts = {
    'A':  'DEGPSTV',
    'R':  'CQGHILKMPSTW',
    'N':  'DHIKSTY',
    'D':  'ANEGHYV',
    'C':  'RGFSWY',
    'E':  'ADQGKV',
    'Q':  'REHLKP',
    'G':  'ARDCESWV',
    'H':  'RNDQLPY',
    'I':  'RNLKMFSTV',
    'L':  'RQHIMFPSWV',
    'K':  'RNEQIMT',
    'M':  'RILKTV',
    'F':  'CILSYV',
    'P':  'ARQHLST',
    'S':  'ARNCGILFPTWY',
    'T':  'ARNIKMPS',
    'W':  'RCGLS',
    'Y':  'NDCHFS',
    'V':  'ADEGILMF'
}

# amino acid scales

# 19 varimax Georgiev scales, each row is a scale, the order of the rows is
# the same as the amino acids in aa_unambiguous_alph
georgiev_scales_mat = numpy.array([
   [ 0.57, -2.8 , -2.02, -2.46,  2.66, -3.08, -2.54,  0.15, -0.39,  3.1 ,
     2.72, -3.89,  1.89,  3.12, -0.58, -1.1 , -0.65,  1.89,  0.79,  2.64],
   [ 3.37,  0.31, -1.92, -0.66, -1.52,  3.45,  1.82, -3.49,  1.  ,
     0.37,  1.88,  1.47,  3.88,  0.68, -4.33, -2.05, -1.6 , -0.09,
    -2.62,  0.03],
   [-3.66,  2.84,  0.04, -0.57, -3.29,  0.05, -0.82, -2.97, -0.63,
     0.26,  1.92,  1.95, -1.57,  2.4 , -0.02, -2.19, -1.39,  4.21,
     4.11, -0.67],
   [ 2.34,  0.25, -0.65,  0.14, -3.77,  0.62, -1.85,  2.06, -3.49,
     1.04,  5.33,  1.17, -3.58, -0.35, -0.21,  1.36,  0.63, -2.77,
    -0.63,  2.34],
   [-1.07,  0.2 ,  1.61,  0.75,  2.96, -0.49,  0.09,  0.7 ,  0.05,
    -0.05,  0.08,  0.53, -2.55, -0.88, -8.31,  1.78,  1.35,  0.72,
     1.89,  0.64],
   [-0.4 , -0.37,  2.08,  0.24, -2.23, -0.  , -0.6 ,  7.47,  0.41,
    -1.18,  0.09,  0.1 ,  2.07,  1.62, -1.82, -3.36, -2.45,  0.86,
    -0.53, -2.01],
   [ 1.23,  3.81,  0.4 , -5.15,  0.44, -5.66,  0.25,  0.41,  1.61,
    -0.21,  0.27,  4.01,  0.84, -0.15, -0.12,  1.39, -0.65, -1.07,
    -1.3 , -0.33],
   [-2.32,  0.98, -2.47, -1.17, -3.49, -0.11,  2.11,  1.62, -0.6 ,
     3.45, -4.06, -0.01,  1.85, -0.41, -1.18, -1.21,  3.43, -1.66,
     1.31,  3.93],
   [-2.01,  2.43, -0.07,  0.73,  2.22,  1.49, -1.92, -0.47,  3.55,
     0.86,  0.43, -0.26, -2.05,  4.2 ,  0.  , -2.83,  0.34, -5.87,
    -0.56, -0.21],
   [ 1.31, -0.99,  7.02,  1.5 , -3.78, -2.26, -1.67, -2.9 ,  1.52,
     1.98, -1.2 , -1.66,  0.78,  0.73, -0.66,  0.39,  0.24, -0.66,
    -0.95,  1.27],
   [-1.14, -4.9 ,  1.32,  1.51,  1.98, -1.62,  0.7 , -0.98, -2.28,
     0.89,  0.67,  5.86,  1.53, -0.56,  0.64, -2.92, -0.53, -2.49,
     1.91,  0.43],
   [ 0.19,  2.09, -2.44,  5.61, -0.43, -3.97, -0.27, -0.62, -3.12,
    -1.67, -0.29, -0.06,  2.44,  3.54, -0.92,  1.27,  1.91, -0.3 ,
    -1.26, -1.71],
   [ 1.66, -3.08,  0.37, -3.85, -1.03,  2.3 , -0.99, -0.11, -1.45,
    -1.02, -2.47,  1.38, -0.26,  5.25, -0.37,  2.86,  2.66, -0.5 ,
     1.57, -2.93],
   [ 4.39,  0.82, -0.89,  1.28,  0.93, -0.06, -1.56,  0.15, -0.77,
    -1.21, -4.79,  1.78, -3.09,  1.73,  0.17, -1.88, -3.07,  1.64,
     0.2 ,  4.22],
   [ 0.18,  1.32,  3.13, -1.98,  1.43, -0.35,  6.22, -0.53, -4.18,
    -1.78,  0.8 , -2.71, -1.39,  2.14,  0.36, -2.42,  0.2 , -0.72,
    -0.76,  1.06],
   [-2.6 ,  0.69,  0.79,  0.05,  1.45,  1.51, -0.18,  0.35, -2.91,
     5.71, -1.43,  1.62, -1.02,  1.1 ,  0.08,  1.75, -2.2 ,  1.75,
    -5.19, -1.31],
   [ 1.49, -2.62, -1.54,  0.9 , -1.15, -2.29,  2.72,  0.3 ,  3.37,
     1.54,  0.63,  0.96, -4.32,  0.68,  0.16, -2.77,  3.73,  2.73,
    -2.56, -1.97],
   [ 0.46, -1.49, -1.71,  1.38, -1.64, -1.47,  4.35,  0.32,  1.87,
     2.11, -0.24, -1.09, -1.34,  1.46, -0.34,  3.36, -5.46, -2.2 ,
     2.87, -1.21],
   [-4.22, -2.57, -0.25, -0.03, -1.05,  0.15,  0.92,  0.05,  2.17,
    -4.18,  1.01,  1.36,  0.09,  2.33,  0.04,  2.67, -0.73,  0.9 ,
    -3.43,  4.77]
])

# 10 BLOSUM62-derived Georgiev scales
georgiev_blosum_scales_mat = numpy.array([
       [ 0.077,  1.014,  1.511,  1.551, -1.084,  1.477,  1.094,  0.849,
         0.716, -1.462, -1.406,  1.135, -0.963, -1.619,  0.883,  0.844,
         0.188, -1.577, -1.142, -1.127],
       [-0.916,  0.189,  0.215,  0.005, -1.112,  0.229,  0.296,  0.174,
         1.548, -1.126, -0.856, -0.039, -0.585,  1.007, -0.675, -0.448,
        -0.733,  2.281,  1.74 , -1.227],
       [ 0.526, -0.86 , -0.046,  0.323,  1.562, -0.67 , -0.871,  1.726,
        -0.802, -0.761, -0.879, -0.802, -0.972, -0.311,  0.382,  0.423,
         0.178,  1.166, -0.582, -0.633],
       [ 0.004, -0.609,  1.009,  0.493,  0.814, -0.355, -0.718,  0.093,
         1.547,  0.382, -0.172, -0.849, -0.528,  0.623, -0.869,  0.317,
        -0.012, -1.61 ,  0.747,  0.064],
       [ 0.24 ,  1.277,  0.12 , -0.991,  1.828, -0.284,  0.5  , -0.548,
         0.35 , -0.599,  0.032,  0.819,  0.236, -0.549, -1.243,  0.2  ,
         0.022,  0.122, -0.119, -0.596],
       [ 0.19 ,  0.195,  0.834,  0.01 , -1.048, -0.075, -0.08 ,  1.186,
        -0.785,  0.276,  0.344,  0.097,  0.365,  0.29 , -2.023,  0.541,
         0.378,  0.239, -0.475,  0.158],
       [ 0.656,  0.661, -0.033, -1.615, -0.742, -1.014, -0.442,  1.213,
         0.655, -0.132,  0.109,  0.213,  0.062, -0.021,  0.845,  0.009,
        -0.304, -0.542,  0.241,  0.014],
       [-0.047,  0.175, -0.57 ,  0.526,  0.379,  0.363,  0.202,  0.874,
        -0.076,  0.198,  0.146,  0.129,  0.208,  0.098, -0.352, -0.797,
        -1.958, -0.398, -0.251,  0.016],
       [ 1.357, -0.219, -1.2  , -0.15 , -0.121,  0.769,  0.384,  0.009,
        -0.186, -0.216, -0.436,  0.176, -0.56 ,  0.433, -0.421,  0.624,
         0.149, -0.349,  0.713,  0.251],
       [ 0.333, -0.52 , -0.139, -0.282, -0.102,  0.298,  0.667,  0.242,
         0.99 ,  0.207, -0.021, -0.85 ,  0.361, -1.288, -0.298, -0.129,
         0.063,  0.499, -0.251,  0.607]
])

def get_georgiev_scale(index, ambiguous=True):
    scale = dict(zip(aa_unambiguous_alph, georgiev_scales_mat[index]))
    return _get_scale(scale, ambiguous)

def get_georgiev_blosum_scale(index, ambiguous=True):
    scale = dict(zip(aa_unambiguous_alph, georgiev_blosum_scales_mat[index]))
    return _get_scale(scale, ambiguous)

def _get_scale(scale, ambiguous):
    if(ambiguous):
        scale.update(_get_non_aa_letter_dict())
    return scale

def _get_non_aa_letter_dict():
    other_letters = aa_ambiguous_alph + aa_special_alph + aa_ter_alph
    return dict(zip(other_letters, len(other_letters) * [0.0]))

georgiev_scales = [get_georgiev_scale(i) for i in xrange(19)]
georgiev_blosum_scales = [get_georgiev_blosum_scale(i) for i in xrange(10)]

# amino acid clusters
# 4 (www.ebi.ac.uk/2can/biology/molecules_small_aatable.html)
# 7 source: taylor85 adjusted version on the url above
# 2 wikipedia aa propensities
aa_subset_dict = {
        'aliphatic_hydrophobic': 'AVLIMPFW',
        'polar_uncharged': 'GSYNQC',
        'acidic': 'ED',
        'basic': 'KRH',
        'aliphatic': 'ILV',
        'aromatic': 'FYWH',
        'charged': 'HKRED',
        'polar': 'YWHKRDETCSNQ',
        'small': 'VCAGTPSDN',
        'tiny': 'AGCST',
        'helix': 'MALEK',
        'sheet': 'YFWTVI'}
aa_subsets = sorted(aa_subset_dict.keys())

# secondary structure
ss_alph = 'CHE'
ss_short = ['col', 'hel', 'str']
ss_name = ['random coil', 'helix', 'strand']

# solvent accessibility
sa_alph = 'BE'
sa_short = ['bur', 'exp']
sa_name = ['buried', 'exposed']

# nucleotides ORDER IS IMPORTANT FOR CODON TABLE
nucleotide_unambiguous_alph = 'TCAG'
codon_aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
ncodon_per_aa = numpy.array([codon_aas.count(l) for l in codon_aas])
nucleotide_ambiguous_alph = 'MRWSYKVHDBN'
nucleotide_alph = nucleotide_unambiguous_alph + nucleotide_ambiguous_alph

# codons
codons_unambiguous = [a+b+c for a in nucleotide_unambiguous_alph
                            for b in nucleotide_unambiguous_alph
                            for c in nucleotide_unambiguous_alph]

codons = [a+b+c for a in nucleotide_alph
                for b in nucleotide_alph
                for c in nucleotide_alph]

map_to_ambiguous_nucleotides = {
    'A':'A', 'C':'C', 'G':'G', 'T':'T', 'U':'T',
    'W':'AT', 'S':'GC', 'M':'AC', 'K':'GT', 'R':'AG', 'Y':'CT',
    'B':'CGT', 'D':'AGT', 'H':'ACT', 'V': 'ACG',
    'N':'ACGT'
}
map_to_unambiguous_aa = dict(zip(list(aa_unambiguous_alph),
                                 list(aa_unambiguous_alph)))
map_to_unambiguous_aa['*'] = '*'
map_to_unambiguous_aa['DN'] = 'B'
map_to_unambiguous_aa['IL'] = 'J'
map_to_unambiguous_aa['GQ'] = 'Z'

def ambiguous_codons(codon):
    return [x + y + z for x in map_to_ambiguous_nucleotides[codon[0]]
                      for y in map_to_ambiguous_nucleotides[codon[1]]
                      for z in map_to_ambiguous_nucleotides[codon[2]]]
def ambiguous_aas(codon):
    return ''.join(sorted(set([codon_table_unambiguous[c]
            for c in ambiguous_codons(codon)])))

def unambiguous_aa(codon):
    return map_to_unambiguous_aa.get(ambiguous_aas(codon), 'X')

codon_table_unambiguous = dict(zip(codons_unambiguous, codon_aas))
codon_table = dict(zip(codons, [unambiguous_aa(c) for c in codons]))

###############################################################################
# Sequence properties
###############################################################################

def hamming_distance(s0, s1):
    assert(len(s0) == len(s1))
    return sum([not s0[i] == s1[i] for i in range(len(s0))])

def dist_one_codons(codon):
    return [c for c in codons_unambiguous if hamming_distance(codon, c) == 1]

def dist_one_amino_acids(codon):
    codons = dist_one_codons(codon)
    return sorted(set([codon_table_unambiguous[c] for c in codons]))

def aa_substitutions():
    tuples = [(l0, l1) for l0 in aa_unambiguous_alph
                       for l1 in aa_unambiguous_alph
                       if not l0 == l1]
    return tuples

def mutations():
    return [(l0, l1) for l0 in nucleotide_unambiguous_alph
                     for l1 in nucleotide_unambiguous_alph
                     if not l0 == l1]

def codon_mutations():
    result = []
    for codon in codons_unambiguous:
        result.extend(codon_mutation(codon))
    return result

def codon_mutation(codon):
    assert(codon in codons_unambiguous)
    return [(codon, c1) for c1 in codons_unambiguous
                        if hamming_distance(codon, c1) == 1]

def single_mutation_aa_substitutions():
    return [(codon_table_unambiguous[c[0]], codon_table_unambiguous[c[1]])
            for c in codon_mutations()]

def possible_single_mutation_aa_substitutions():
    return sorted(set([s for s in single_mutation_aa_substitutions()
                  if not(s[0] == s[1] or (s[0] == '*' or s[1] == '*'))]))

def impossible_single_mutation_aa_substitutions():
    return sorted(set(aa_substitutions()) -
                  set(possible_single_mutation_aa_substitutions()))

def single_mutation_aa_substitution_stats():
    subs = single_mutation_aa_substitutions()
    no_sub = 0
    stop_sub = 0
    aa_subs = []
    for s in subs:
        if(s[0] == s[1]):
            no_sub += 1
        elif(s[0] == '*' or s[1] == '*'):
            stop_sub += 1
        else:
            aa_subs.append(s)
    print ''
    print 'No substitution: %i' % (no_sub)
    print 'Stop codon substitution: %i' % (stop_sub)
    print 'Amino acid substitution: %i' % (len(aa_subs))
    print 'TOTAL: %i' % (len(subs))
    print ''
    possible_subs = set(aa_subs)
    impossible_subs = set(aa_substitutions()) - possible_subs
    print 'Possible amino acid substitutions: %i' % (len(possible_subs))
    print 'Impossible amino acid substitutions: %i' % (len(impossible_subs))
    print 'TOTAL: %i' % (len(aa_substitutions()))
    #print impossible_subs
    print ''

###############################################################################
# Sequence operations (not yet optimized for speed...)
###############################################################################


def translate(orf):
    '''
    This function translates a possibly ambiguous ORF nucleotide sequence into
    an amino acid protein sequence.

    Args:
        orf (str): The open reading frame (nucleotide) sequence.

    Returns:
        str The translated protein (amino acid) sequence.

    Raises:
        ValueError: if the orf length is not a multiple of 3.
        ValueError: if orf is not an (possibly ambiguous) nucleotide sequence.

    Translation of a the start of some random protein can be done with:

    >>> translate('ATGTTTAGTAACAGACTACCACCTCCAAAA')
    'MFSNRLPPPK'
    '''
    if not(len(orf) % 3 == 0):
        raise ValueError('ORF sequence length is not a multiple of 3.')
    if not(is_nucleotide_sequence(orf)):
        raise ValueError('ORF sequence is not a nucleotide sequence.')

    return ''.join([codon_table[orf[i:i+3]] for i in xrange(0, len(orf), 3)])

def seq_count(seq, alph):
    '''
    This function counts letter occurance in seq for each letter in alph.

    Args:
        seq (str): The sequence of which the letters will be counted.
        alph (str): The letters that will be counted in seq

    Returns:
        numpy.array List with letter counts in the order of alph.

    >>> seq_count('AABBCBBACB', 'ABC')
    array([3, 5, 2])
    >>> seq_count('', 'ABC')
    array([0, 0, 0])
    >>> seq_count('ABC', '')
    array([], dtype=int64)
    '''
    return numpy.array([seq.count(l) for l in alph], dtype=int)

def seq_composition(seq, alph):
    '''
    '''
    return seq_count(seq, alph) / float(len(seq))

def aa_count(protein):
    return seq_count(protein, aa_unambiguous_alph)

def aa_composition(protein):
    return seq_composition(protein, aa_unambiguous_alph)

def ss_composition(protein_ss):
    return seq_composition(protein_ss, ss_alph)

def sa_composition(protein_sa):
    return seq_composition(protein_sa, sa_alph)

def state_subseq(seq, state_seq, state_letter):
    assert(len(seq) == len(state_seq))
    return ''.join([l if state_seq[i] == state_letter else ''
                      for i, l in enumerate(seq)])

def state_subseq_composition(seq, state, seq_alph, state_alph):
    result = []
    for l in state_alph:
        result.extend(seq_composition(state_subseq(seq, state, l), seq_alph))
    return result

def ss_aa_composition(protein, ss):
    return state_subseq_composition(protein, ss, aa_unambiguous_alph, ss_alph)

def sa_aa_composition(protein, sa):
    return state_subseq_composition(protein, sa, aa_unambiguous_alph, sa_alph)

def aa_cluster_count(protein):
    counts = dict(zip(aa_unambiguous_alph, aa_count(protein)))
    return numpy.array([sum([comp[l] for l in aa_subset_dict[subset]])
            for subset in aa_subsets])

def aa_cluster_composition(protein):
    comp = dict(zip(aa_unambiguous_alph, aa_composition(protein)))
    return numpy.array([sum([comp[l] for l in aa_subset_dict[subset]])
            for subset in aa_subsets])

def window_seq(seq, window_size, overlapping=False):
    '''
    If the length of the sequence is not a multiple of the window size, the
    last letters are NOT returned.
    >>> from util import sequtil
    >>> s = 'AAAACCACCAAAA'
    >>> sequtil.window_seq(s, 3)
    ['AAA', 'ACC', 'ACC', 'AAA']
    '''
    if(window_size < 2):
        return list(seq)
    else:
        start = 0
        stop = len(seq) - window_size + 1
        step = window_size
        if(overlapping):
            step = 1
        return [seq[i:i + window_size] for i in range(start, stop, step)]

filter_cache = {}
def convolution_filter(window=9, edge=0.0):  # type not implemented yet
    '''
    Returns triangular convolution filter. The filter values add up to 1.0.
    '''

    if((window, edge) in filter_cache.keys()):
        return filter_cache[(window, edge)]

    if(window % 2 == 0):
        raise ValueError('Window must be an uneven number.')
    if(window < 3):
        raise ValueError('Window must be 3 or larger.')
    if(edge < 0.0 or edge > 1.0):
        raise ValueError('The edge parameter must be in the range 0.0 to 1.0.')

    if(edge == 1.0):
        result = numpy.ones(window) / window
        filter_cache[(window, edge)] = result
        return result
    else:
        result = numpy.ones(window)
        num = window / 2
        #diff = 1.0 - edge
        #step = diff / distance
        #forw = numpy.arange(edge, 1.0, step)
        #backw = forw[::-1]
        forw = numpy.linspace(edge, 1.0, num, endpoint=False)
        result[:num] = forw
        result[-num:] = forw[::-1]
        #forw = numpy.append(forw, 1.0)
        #filt = numpy.append(forw, backw)
        result = result / result.sum()
        filter_cache[(window, edge)] = result
        return result

def seq_signal_raw(protein, scale):
    return [scale[letter] for letter in protein]

def seq_signal(protein, scale, window=9, edge=0.0):

    if(window > len(protein) or window < 1):
        raise ValueError('1 <= window <= length protein.')

    # obtain raw signal
    signal = seq_signal_raw(protein, scale)

    # return the raw signal if window size is one
    if(window == 1):
        return signal

    # otherwise return convolved signal
    else:
        conv = convolution_filter(window, edge)
        return numpy.convolve(signal, conv, 'valid')

def codon_count(orf):
    wseq = window_seq(orf, 3, overlapping=False)
    return numpy.array([wseq.count(c) for c in codons_unambiguous])

def codon_composition(orf):
    return codon_count(orf) / float(len(orf))

def codon_usage(orf):
    aa_c_dict = dict(zip(aa_unambiguous_alph, aa_count(translate(orf))))
    # change 0 to 1, to prevent / 0 (answer will still always be 0 (0/1))
    aa_c = numpy.array([float(aa_c_dict.get(a, 0))
                           if aa_c_dict.get(a, 0) > 0 else 1.0
                           for a in codon_aas])
    codon_c = codon_count(orf)
    return codon_c / aa_c

###############################################################################
# Sequence checkers
###############################################################################
def is_amino_acid_sequence(sequence):
    if(set(sequence).issubset(set(aa_alph))):
        return True
    else:
        print sorted(set(sequence))
        print sorted(set(aa_alph))
        return False

def is_unambiguous_amino_acid_sequence(sequence):
    if(set(sequence).issubset(set(aa_unambiguous_alph))):
        return True
    return False

def is_nucleotide_sequence(sequence):
    if(set(sequence).issubset(set(nucleotide_alph))):
        return True
    return False

def is_unambiguous_nucleotide_sequence(sequence):
    if(set(sequence).issubset(set(nucleotide_unambiguous_alph))):
        return True
    return False

def is_sec_struct_sequence(sequence):
    if(set(sequence).issubset(set(ss_alph))):
        return True
    return False

def is_solv_access_sequence(sequence):
    if(set(sequence).issubset(set(sa_alph))):
        return True
    return False
