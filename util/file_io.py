'''
Created on Sep 10, 2011

@author: Bastiaan van den Berg
'''

'''
class Fasta(object):

    def __init__(self, fasta_f):

    def __iter__(self):
        pass

    def next(self):
        pass
'''

# This is a generator function
def read_fasta(f, filter_ids=None):
    '''
    '''

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'r')

    # initialize sequence id and string to an empty string
    seq_id = ""
    seq_str = ""

    # iterate over each line in the fasta file
    for line in handle:
        
        if(seq_id == "" and seq_str == ""):
            if(line[0] == ">"):
                seq_id = line.split()[0][1:]
            elif(line[0] == '#'):
                pass # comment
            elif(line.strip()):
                # non-empty line...
                print( line.strip())
                raise(Exception, "Error in fasta file")
        else:
            if(line == "" or line[0] == ">"):
                if(filter_ids == None or seq_id in filter_ids):
                    yield (seq_id, seq_str)
                seq_str = ""
                if(line[0] == ">"):
                    seq_id = line.split()[0][1:]
                else:
                    seq_id = ""
            else:
                seq_str += line.strip()

    # return the last sequence (not if the file was empty)
    if not(seq_id == ""):
        if(filter_ids == None or seq_id in filter_ids):
            yield (seq_id, seq_str)

    # close file if we opened it
    if not(type(f) == file):
        handle.close()

def write_fasta(f, seqs):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'w')

    for s in seqs:
        handle.write('>' + s[0] + '\n')
        for i in range(0, len(s[1]), 70):
            handle.write(s[1][i:i + 70] + '\n')
        handle.write('\n')
        handle.flush()

    # close file if we opened it
    if not(type(f) == file):
        handle.close()

def write_mutation(f, mutations):
    write_tuple_list(f, mutations)
    
def read_mutation(f):
    return read_tuple_list(f, (str, int, str, str))

def read_labeling(f):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'r')

    class_names = handle.readline().split()

    label_dict = {}
    for line in handle:
        tokens = line.split()
        label = int(tokens[1])
        assert(label < len(class_names))
        label_dict[tokens[0]] = label

    if not(type(f) == file):
        handle.close()

    return (label_dict, class_names)

def write_labeling(f, object_ids, labels, class_names):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'w')

    handle.write('%s\n' % ('\t'.join(class_names)))

    for (obj, lab) in zip(object_ids, labels):
        handle.write('%s\t%s\n' % (obj, lab))

    # close file if we opened it
    if not(type(f) == file):
        handle.close()

def read_propka30(filename):

    # values to be read    
    feph = [] # 15 times free energy per pH
    chphf = [] # 15 charge per pH folded
    chphu = [] # 15 charge per pH unfolded
    femin = 100000.0
    feminph = -1.0
    pif = -1.0
    piu = -1.0

    # parse status
    first = True
    fe = False
    fe_count = 0
    ch = False
    ch_count = 0
    max_count = 15

    with open(filename, 'r') as fin:

        for line in fin:
            tokens = line.split()
            if(first):
                assert(tokens[0] == 'propka3.0,' and tokens[2] == '182')
                first = False
            if(tokens): 
                if(fe and fe_count < max_count):
                    feph.append(float(tokens[1]))
                    fe_count += 1
                if(ch and ch_count < max_count):
                    chphf.append(float(tokens[1]))
                    chphu.append(float(tokens[2]))
                    ch_count += 1
                if(tokens[0] == 'Free'):
                    fe = True
                elif(tokens[0] == 'pH'):
                    ch = True
                if(tokens[0] == 'The' and tokens[1] == 'pI'):
                    pif = float(tokens[3])
                    piu = float(tokens[6])
                if(tokens[0] == 'The' and tokens[1] == 'pH'):
                    feminph = float(tokens[6])
                    femin = float(tokens[13])

    assert(len(feph) == max_count)
    assert(len(chphf) == max_count)
    assert(len(chphu) == max_count)
    assert(not femin == 100000.0)
    assert(not feminph == -1.0)
    assert(not pif == -1.0)
    assert(not piu == -1.0)

    # bit tricky... return as dict???
    return((feph, chphf, chphu, femin, feminph, pif, piu))

# This is a generator function
def read_ids(f):
    
    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'r')
    
    for line in handle:
        tokens = line.split()
        yield(tokens[0])

    # close file if we opened it
    if not(type(f) == file):
        handle.close()

def write_ids(f, ids):
    
    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'w')
    
    for uid in ids:
        handle.write('%s\n' % (uid))
    
    # close file if we opened it
    if not(type(f) == file):
        handle.close()

def read_dict(handle, value_type, num_cols=1):
    result = {}
    for line in handle:
        tokens = line.split()
        if(num_cols > 1):
            end = num_cols + 2
            result[tokens[0]] = tuple([value_type(i) for i in tokens[1:end]])
        else:
            result[tokens[0]] = value_type(tokens[1])
    return result

def write_dict(handle, d):
    for key in d.keys():
        handle.write('%s\t%s\n' % (key, str(d[key])))

# use eval instead of passing list of types???
def read_tuple_list(f, types):
    
    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'r')
    
    tuples = []
    for line in handle:
        tokens = line.split()
        row = []
        for index in xrange(len(types)):
            row.append(types[index](tokens[index]))
        tuples.append(tuple(row))
    
    # close file if we opened it
    if not(type(f) == file):
        handle.close()
    
    return tuples

def write_tuple_list(f, tuple_list):
    
    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'w')
    
    for tup in tuple_list:
        for item in tup:
            handle.write('%s\t' % (str(item)))
        handle.write('\n')
    
    # close file if we opened it
    if not(type(f) == file):
        handle.close()
