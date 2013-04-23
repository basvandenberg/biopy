import itertools
import numpy
import random
from operator import mul

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
aa_ambiguous_name = ['aspartic acid or asparagine', 'leucine or isoleucine', 
                      'unknown amino acid', 'glutamic acid or glutamine']

aa_special_alph = 'UO'
aa_special_short = ['sec', 'pyl']
aa_special_name = ['selenocysteine', 'pyrralysine']

aa_ter_alph = '*'
aa_ter_short = ['ter']
aa_ter_name = ['terminal']

aa_alph = aa_unambiguous_alph + aa_ambiguous_alph +\
          aa_special_alph + aa_ter_alph
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
       [ 0.57, -2.8 , -2.02, -2.46,  2.66, -3.08, -2.54,  0.15, -0.39,
         3.1 ,  2.72, -3.89,  1.89,  3.12, -0.58, -1.1 , -0.65,  1.89,
         0.79,  2.64],
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

################################################################################
# Sequence properties
################################################################################

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

################################################################################
# Sequence operations (not yet optimized for speed...)
################################################################################

# translate (ambiguous) orf into (ambiguous) amino acid sequence
def translate(orf):
    return ''.join([codon_table[orf[i:i+3]] for i in xrange(0, len(orf), 3)])

def seq_count(seq, alph):
    return numpy.array([seq.count(l) for l in alph])

def seq_composition(seq, alph):
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
    If the length of the sequence is not a multiple of the window size, the last
    letters are NOT returned.
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






class Seq(str):
    """
    classdocs

    A sequence may only contain letters from the alphabet. The sequence does
    not have to contain all the letters from the alphabet.
    """
    
    def __new__(cls, value, *args, **keywargs):
        return str.__new__(cls, value)

    def __init__(self, value, seq_id, alph, ambiguous=None):

        self.seq_id = seq_id
        self.alph = alph

        if(ambiguous):
            self.alph_unamb = alph
            self.alph = self.alph + ambiguous

        self.letter_count = None

        if not(set(value) <= set(self.alph)):
            raise ValueError('Incorrect letters in sequence.')

    def shuffled(self, window_size=1):
        '''
        Returns a sequence with the same content, but shuffled.

        >>> from bio import seq
        >>> s = seq.Seq('AAAACCACCAAA', 'AC')
        >>> s_shuf = s.shuffled()
        >>> len(s_shuf)
        12
        >>> s_shuf.count('A')
        8
        >>> s_shuf.count('C')
        4
        >>> s_shuf = s.shuffled(window_size=3)
        >>> len(s_shuf)
        12
        >>> s_shuf.window_seq(window_size=3).count('AAA')
        2
        >>> s_shuf.window_seq(window_size=3).count('ACC')
        2
        '''
        l = self.window_seq(window_size)
        random.shuffle(l)
        return Seq(''.join(l), self.alphabet)

    def int_seq(self, window_size=1, nalphabet=None):
        return [self.int_map(window_size)[letter] for letter in self.window_seq(window_size)]

    def from_int_seq(self, int_seq, window_size=1):
        return [self.rev_int_map(window_size)[i] for i in int_seq]

    def rev_int_map(self, window_size=1):
        imap = self.int_map(window_size=window_size)
        return dict((v,k) for k, v in imap.iteritems())

    def int_map(self, window_size=1):
        '''
        >>> from bio import seq
        >>> s = seq.Seq('', 'AC')
        >>> s.int_map()
        {'A': 0, 'C': 1}
        >>> s.int_map(window_size=2)
        {'AA': 0, 'CC': 3, 'AC': 1, 'CA': 2}
        '''
        mapping = {}
        letters = sorted(self.alph_product(window_size))
        for i in range(len(letters)):
            mapping[letters[i]] = i
        return mapping

    def alph_product(self, window_size):
        return [''.join(i) for i in itertools.product(self.alphabet, repeat=window_size)]
    
    def window_seq(self, window_size, overlapping=False):
        '''
        If the length of the sequence is not a multiple of the window size, the last
        letters are NOT returned.
        >>> from bio import seq
        >>> s = seq.Seq('AAAACCACCAAAA', 'AC')
        >>> s.window_seq(3)
        ['AAA', 'ACC', 'ACC', 'AAA']
        '''
        if(window_size < 2):
            return list(self)
        else:
            start = 0
            stop = len(self) - window_size + 1
            step = window_size
            if(overlapping):
                step = 1
            return [self[i:i + window_size] for i in range(start, stop, step)]

    def mapped_seq(self, mapping, window_size=1):
        '''
        '''
        mapseq = ''.join([mapping[letter] for letter in self.window_seq(window_size)])
        alph = ''.join(set(''.join(mapping.values())))
        return Seq(mapseq, alph)

    
    
    ### Feature calculation related, not used for a while...

    def count_letters(self):
        """ Returns letter count of each letter in the alphabet.

        This method returns a dictionary that maps each letter in the alphabet
        to the number of times it occurs in this sequence.

        >>> from bio import seq
        >>> s = seq.Seq('ueeueeuaue', 'aoeu')
        >>> letter_count = s.count_letters()
        >>> letter_count['a']
        1
        >>> letter_count['o']
        0
        >>> letter_count['e']
        5
        >>> letter_count['u']
        4
        """
        # check if letters are already counted
        if(self.letter_count == None):
            # create dict with each letter in ALPHABET initialized to zero
            self.letter_count = dict(zip(self.alphabet,
                                         len(self.alphabet) * [0]))
            # count occurances of each letter in seq_str
            for letter in self:
                self.letter_count[letter] += 1
        # return letter count dictionary
        return self.letter_count
    
    def cl(self, letter):
        return self.count_letters()[letter]
    
    def count_letter_set(self, letters):
        return sum(map(self.cl, letters))
    
    def composition(self, letters):
        """ Returns the composition for the given letters.

        pre: letters must be a subset of self.alphabet.

        >>> from bio import seq
        >>> s = seq.Seq('ueeueeuuue', 'aoeu')
        >>> s.composition(s.alphabet)
        {'a': 0.0, 'u': 0.5, 'e': 0.5, 'o': 0.0}
        >>> s.composition('ae')
        {'a': 0.0, 'e': 0.5}
        >>> s.composition('u')
        0.5
        """
        if(len(letters) == 1):
            return self._letter_composition(letters[0])
        else:
            return dict(zip(letters, map(self._letter_composition, letters)))

    def _letter_composition(self, letter):
        if(len(self) > 0):
            return float(self.count_letters()[letter]) / len(self)
        else:
            return 0.0
    
    # maybe remove and use composition of MappedSeq(Seq) class instead.
    def composition_set(self, letters):
        """ Returns the combined composition of a set of letters.

        >>> from bio import seq
        >>> s = seq.Seq('ueeueeuaae', 'aoeu')
        >>> s.composition_set('au')
        0.5
        """
        return sum(self.composition(letters).values())

    def signal(self, scale, filt_win=1):
        """ Returns a protein signal for the give amino acid scale. The signal
        will be filtered if a filter is supplied.

        >>> from bio import seq
        >>> scale = {'a': 0.25, 'o': 0.5, 'e': 0.75, 'u': 0.0}
        >>> s = seq.Seq('ueeueeuaue', 'aoeu')
        >>> sig = s.signal(scale)
        >>> sig.signal_raw
        [0.0, 0.75, 0.75, 0.0, 0.75, 0.75, 0.0, 0.25, 0.0, 0.75]
        """
        return SeqSignal(self, scale, filt_win)

    def encode(self, encoding):
        '''
        similar te signal... used for weighted degree string kernel
        '''
        result = []
        [result.extend(encoding[l]) for l in self]
        return result
        
    def pair_occ(self, pairs, gap):
        """
        >>> from bio import seq
        >>> s = seq.Seq('aoaoaoaoaoao', 'ao')
        >>> s.pair_occ('aa', 0)
        -13.815510557964274
        >>> s.pair_occ(['oo', 'ao'], 0)
        {'oo': -13.815510557964274, 'ao': 0.78015855754957497}
        """
        if not(type(pairs) == list):
            return self._pair_occ(pairs, gap)
        else:
            result = {}
            for pair in pairs:
                result[pair] = self._pair_occ(pair, gap)
            return result

    def _pair_occ(self, pair, gap):
        di_comp = self.composition_pair(pair, gap)
        di_prob = self.pair_prob(pair)
        frac = 0.000001
        if not(di_comp == 0):
            frac = di_comp / di_prob
        return numpy.log(frac)

    def composition_pair(self, pairs, gap):
        """
        >>> from bio import seq
        >>> s = seq.Seq('aoaoaoaoaoa', 'ao')
        >>> s.composition_pair(['aa', 'oo', 'oa', 'ao'], 0)
        {'aa': 0.0, 'oo': 0.0, 'ao': 0.5, 'oa': 0.5}
        >>> s.composition_pair(['aa', 'ao'], 0)
        {'aa': 0.0, 'ao': 0.5}
        >>> s.composition_pair('aa', 0)
        0.0
        """
        if not(type(pairs) == list):
            return self._pair_composition(pairs, gap)
        else:
            result = {}
            for pair in pairs:
                result[pair] = self._pair_composition(pair, gap)
            return result
    
    def count_pairs(self, gap):
        """
        >>> from bio import seq
        >>> s = seq.Seq('ueeueeuuue', 'aoeu')
        >>> pdict = s.count_pairs(1)
        >>> pdict['ue']
        3
        >>> pdict['uu']
        1
        """
        pair_count = {}
        for i in range(self.num_pairs(gap)):
            pair = self[i] + self[i + gap + 1]
            pair_count[pair] = pair_count.setdefault(pair, 0) + 1
        return pair_count

    def _pair_composition(self, pair, gap):
        try:
            return float(self.count_pairs(gap)[pair]) / self.num_pairs(gap)
        except (KeyError, ZeroDivisionError):
            return 0.0
    
    def pair_prob(self, pair):
        """
        >>> from bio import seq
        >>> s = seq.Seq('aoaoaoaoaoao', 'ao')
        >>> s.pair_prob('ao')
        0.25
        >>> s.pair_prob('aa')
        0.25
        """
        if(pair[0] == pair[1]):
            return self.composition(pair[0]) ** 2
        else:
            return reduce(mul, self.composition(pair).values())

    def num_pairs(self, gap):
        """
        >>> from bio import seq
        >>> s = seq.Seq('aoeuaoeu', 'aoeu')
        >>> s.num_pairs(0)
        7
        >>> s.num_pairs(4)
        3
        """
        return self.num_words(2 + gap)

    def count_triplets(self, gap1, gap2):
        """
        >>> from bio import seq
        >>> s = seq.Seq('ueeueeuuue', 'aoeu')
        >>> pdict = s.count_triplets(0,1)
        >>> pdict['ueu']
        2
        >>> pdict['eue']
        1
        """
        triplet_count = {}
        for i in range(self.num_triplets(gap1, gap2)):
            triplet = self[i] + \
                      self[i + gap1 + 1] +\
                      self[i + gap1 + gap2 + 2]
            triplet_count[triplet] = triplet_count.setdefault(triplet, 0) + 1
        return triplet_count

    def composition_triplet(self, triplets, gap1, gap2):
        """
        >>> from bio import seq
        >>> s = seq.Seq('aaaaaaaooooo', 'ao')
        >>> tcomp = s.composition_triplet(['aaa', 'aoa'], 0, 0)
        >>> tcomp['aaa']
        0.5
        >>> tcomp['aoa']
        0.0
        >>> s.composition_triplet('aaa', 0, 0)
        0.5
        """
        gap1 = int(gap1)
        gap2 = int(gap2)
        if(type(triplets) == str):
            return self._triplet_composition(triplets, gap1, gap2)
        else:
            result = {}
            for tri in triplets:
                result[tri] = self._triplet_composition(tri, gap1, gap2)
            return result

    def _triplet_composition(self, triplet, gap1, gap2):
        try:
            return float(self.count_triplets(gap1, gap2)[triplet]) /\
                         self.num_triplets(gap1, gap2)
        except (KeyError, ZeroDivisionError):
            return 0.0

    def composition_triplet_set(self, triplets, gap1, gap2):
        """
        >>> from bio import seq
        >>> s = seq.Seq('aaaaaaaooooo', 'ao')
        >>> s.composition_triplet_set(['aaa', 'aoa'], 0, 0)
        0.5
        """
        return sum(self.composition_triplet(triplets, gap1, gap2).values())
    
    def num_triplets(self, gap1, gap2):
        """
        >>> from bio import seq
        >>> s = seq.Seq('aoeuaoeu', 'aoeu')
        >>> s.num_triplets(0, 0)
        6
        >>> s.num_triplets(1, 2)
        3
        """
        return self.num_words(3 + gap1 + gap2)

    def num_words(self, word_length, window=1):
        return len(self) / window - word_length + 1

class StateSeq(Seq):
    """
    classdocs
    """

    def __init__(self, value, alphabet, parent):
        super(StateSeq, self).__init__(value, alphabet)
        self.parent = parent
        self._check_state_seq()

    def _check_state_seq(self):
        if not(len(self) == len(self.parent)):
            raise(Exception, "State and parent seq not same length.")

    def parent_int_seq(self):
        self_seq = self.int_seq()
        parent_seq = self.parent.int_seq()
        parent_alph_size = len(self.parent.alphabet)
        return [parent_seq[i] + parent_alph_size * self_seq[i] for i in range(len(self))]

    def mapped_parent_int_seq(self, mapping):
        self_seq = self.int_seq()
        parent_mapped_seq = self.parent.mapped_seq(mapping)
        parent_mapped_seq_ints = parent_mapped_seq.int_seq()
        parent_alph_size = len(parent_mapped_seq.alphabet)
        return [parent_mapped_seq_ints[i] + parent_alph_size * self_seq[i] for i in range(len(self))]
    
    def parent_seq(self, state_letter):
        """
        >>> from bio import seq
        >>> s = seq.Seq('auuueeueoo', 'aoeu')
        >>> r = seq.StateSeq('hhhhtttthh', 'ht', s)
        >>> r.parent_seq('h')
        'auuuoo'
        >>> r.parent_seq('t')
        'eeue'

        """
        return Seq("".join([self.parent[i] for i in range(len(self))
                            if self[i] == state_letter]),
                            self.parent.alphabet)

    def parent_seqs(self, state_letters, int_seqs=True):
        """
        >>> from bio import seq
        >>> s = seq.Seq('auuueeueoo', 'aoeu')
        >>> r = seq.StateSeq('hhhhtttthh', 'ht', s)
        >>> r.parent_seqs('h')
        ['auuu', 'oo']
        >>> r.parent_seqs('t')
        ['eeue']
        """
        
        if(int_seqs):
            parent_seq = self.parent.int_seq()
        else:
            parent_seq = self.parent
                
        seqs = []
        current_seq = []
        for i in range(len(self)):
            if(self[i] in state_letters):
                current_seq.append(parent_seq[i])
            elif not(current_seq == []):
                seqs.append(current_seq)
                current_seq = []
        if not(current_seq == []):
            seqs.append(current_seq)

        if not(int_seqs):
            seq_copy = []
            for s in seqs:
                seq_copy.append(''.join(s))
            seqs = seq_copy

        return seqs


class SeqSignal(object):

    def __init__(self, seq, scale, filt_win=1, scalewin_size=1):
        '''
        Constructor
        '''
        self.seq = seq
        self.scale = scale
        self.filt_win = filt_win
        self.scalewin_size = scalewin_size
        self.signal_raw = self._signal_raw()
        self.signal_smoothed = self._signal_smoothed()

    def _signal_raw(self):
        return map(lambda letter: self.scale[letter],
                   self.seq.window_seq(self.scalewin_size, False))

    def _signal_smoothed(self):
        if(self.filt_win < 3):
            return self.signal_raw
        else:
            return numpy.convolve(self.signal_raw,
                                  self.convolution_filter(self.filt_win),
                                  'valid')

    def smooth(self, convolution_filter):
        """

        """
        self.convolution_filter = convolution_filter
        self.signal_smoothed = self._signal_smoothed()
        return self.smoothed()

    def avg(self):
        """

        >>> from bio import seq
        >>> scale = {'a': 0.25, 'o': 0.5, 'e': 0.75, 'u': 0.0}
        >>> s = seq.Seq('uuaaaaoo', 'aoeu')
        >>> s.signal(scale).avg()
        0.25
        """
        if(len(self.seq) == 0):
            return 0.0
        else:
            return float(sum(self.signal_smoothed)) / len(self.seq)

    def auc(self, threshold):
        """

        >>> from bio import seq
        >>> scale = {'a': 0.25, 'o': 0.5, 'e': 0.75, 'u': 0.0}
        >>> s = seq.Seq('uuaaaauoo', 'aoeu')
        >>> (above, below) = s.signal(scale).auc(0.25)
        >>> above
        0.5
        >>> below
        0.75
        """
        area_above = 0.0
        area_below = 0.0
        for value in self.signal_smoothed:
            if value > threshold:
                area_above += value - threshold
            elif value < threshold:
                area_below += threshold - value
            else:
                pass
        return (area_above, area_below)

    def norm_auc(self, threshold, location):
        if(len(self.seq) > 0):
            if(location == 'top'):
                return float(self.auc(threshold)[0]) / len(self.seq)
            elif(location == 'bot'):
                return float(self.auc(threshold)[1]) / len(self.seq)
            else:
                raise ValueError('Location should be\'top\' or\'bot\'.')
        else:
            return 0.0

#    def plot(self, handle=None, c='#2465a4'):
#        start = 0
#        if not(self.convolution_filter is None):
#            start = len(self.convolution_filter) / 2
#        y = self.signal_smoothed
#        x = numpy.arange(start, start + len(y), 1)
#        pyplot.plot(x, y, color=c)
#        pyplot.xlim(0, 1260)
#        pyplot.ylim(0.0, 1.0)
#        if(handle is None):
#            pyplot.show()
#        else:
#            pyplot.savefig(handle, dpi=600)

