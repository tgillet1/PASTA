'''
Provides the ability to model and represent various concrete objects that
are to be used throughout application runtime.
'''

from collections import Counter

default_nodetypes = {'A':'A','C':'C','T':'T'}

# Reads in a file containing node types (A,C,T) and a string of specific characters of that type
# Format for a given line should be: "<node-type>:<character string>", for example: "C:BRPD", or the simple case: "C:C"
def parse_node_types(fname):
    node_types = {}
    for line in open(fname):
        line = line.strip()
        if len(line) == 0:
            break
        else:
            vals = line.split(':')
            node_types[vals[0]] = vals[1]
    return node_types

def generate_identity_matrix(nodetypes=default_nodetypes):
    submat = {}
    ac_types = ['A','C']
    for char1 in nodetypes.keys():
        nodetype1 = nodetypes[char1]
        for char2 in nodetypes.keys():
            nodetype2 = nodetypes[char2]
            if nodetype1 == nodetype2 or nodetype1 in ac_types and nodetype2 in ac_types:
                submat[(char1,char2)] = 1
                submat[(char2,char1)] = 1
            else:
                submat[(char1,char2)] = -40
                submat[(char2,char1)] = -40
    return submat


class NeuriteSequence():
    ''' 
    A NeuriteSequence object is simply a FASTA object but solely references
    neuronal sequences.
    @param sequence: input sequence string.
    @param header: header for the respective sequence.
    '''
    def __init__(self, name, seq):
        self.seq = seq
        self.name = name
        
    def get_length(self):
        ''' 
        Returns the length of the neurite-sequence object.
        @return: integer referencing sequence length.
        '''
        return len(self.seq)
    
    def get_name(self):
        '''
        Returns the name of the neurite-sequence object.
        @return: string referencing sequence object.
        '''
        return self.name
    
    def get_sequence(self):
        ''' 
        Returns the current sequence as a string.
        @return: string referencing the current neuronal sequence.
        '''
        return self.seq
    
    def __str__(self):
        return '{ ' + self.get_name() + '; ' + self.get_sequence() +\
            '; ' + str(self.get_length()) + ' bases }'
    
    def __repr__(self):
        return self.__str__()

tree_sequence_types = ['complete_tree','incomplete_tree','multiple_trees_complete','multiple_trees_incomplete']

def tree_sequence_type(seq,node_types=default_nodetypes):
    a_nodes = 0
    #AStack.append(-1) # Begin the AStack with a sentinal value
    complete_trees = 0
    tree_complete = True
    # Visit each character in the sequence
    for index in range(0,len(seq)):
        tree_complete = False
        # If the current node is an A-node
        if seq[index] in node_types['A']:
            # push it onto the top of the stack
            a_nodes += 1
 
        # If the current node is a T-node
        elif seq[index] in node_types['T']:
            if a_nodes == 0:
                complete_trees += 1
                tree_complete = True
            else:
                a_nodes -= 1

    if complete_trees == 0:
        return 'incomplete_tree'
    elif tree_complete:
        if complete_trees == 1:
            return 'complete_tree'
        else: # more than one complete tree
            return 'multiple_trees_complete'
    else:
        return 'multiple_trees_incomplete'

    
class ConsensusSequence(NeuriteSequence):
    ''' 
    Encapsulates a consensus sequence so it can be parsed and analyzed using 
    domain-analysis mode.
    '''
#    def __init__(self, name, seq, threshold=None, conservation_score=None):
#        super(ConsensusSequence, self).__init__(name, seq)
#        self.num_gaps = 0
#        self.threshold = threshold
#        self.conservation_score = conservation_score
#        self.ungapped_consensus = seq.replace('-','')
 
    def __init__(self,consensus,threshold_ratio,height,threshold_type='percent',average_conservation=None,name=''):
        super(ConsensusSequence, self).__init__(name,consensus)
        self.raw_consensus = self.seq
        self.consensus = self.seq.replace('-','')
        self.height = height
        self.threshold_ratio = threshold_ratio
        self.threshold_type = threshold_type
        self.average_conservation = average_conservation

    def get_consensus_length(self):
        return len(self.consensus)

    def get_threshold_count(self):
        if self.threshold_type == 'sqrt':
            threshold_count = self.threshold_ratio*math.sqrt(self.height) # consensus threshold                                                        
        else: # Default to percent                                                                                                                     
            threshold_count = self.threshold_ratio*self.height # consensus threshold                                                                   
        return threshold_count

    def get_threshold_percent(self):
        if self.threshold_type == 'sqrt':
            threshold_percent = self.threshold_ratio*math.sqrt(self.height)/self.height
        else: # Default to percent                                                                                                                     
            threshold_percent = self.threshold_ratio
        return threshold_percent

def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if not length % 2:
        return (sorts[int(length / 2)] + sorts[int(length / 2 - 1)]) / 2.0
    return sorts[int(length / 2)]
       
class MultipleSequenceAlignment():
    '''
    An object class containing a complete multiple sequence alignment, including component
    sequences with MSA gaps, composite sequence, and optionally one or more consensus
    sequences. Also contains methods for producing additional consensuses and MSA statistics.
    '''

    def __init__(self,composite,alignments={},node_types=default_nodetypes,consensuses={}):
        self.composite = composite
        self._alns = alignments
        self.node_types = node_types
        self.consensuses = consensuses

        self.composite_node_types = NeuriteSequence('PWM Node Types',''.join(list([self.node_types[char] for char in self.composite])))

        if len(alignments) > 0:
            self._build_stats()
        else:
            # Flag noting whether col_counts needs to be updated by calling _build_stats                                              
            self._reset_stats = 1

    def set_alignments(self,alignments):
        self._alns = alignments
        self._reset_stats = 1

    def add_alignment(self,name,alignment):
        self._alns[name] = alignment
        self._reset_stats = 1

    def get_alignment(self,name):
        return self._alns[name]

    def get_alignments(self):
        return self._alns.copy()

    def add_consensus(self,threshold,consensus,threshold_type='percent'):
        self.consensuses[str(threshold)+threshold_type] = consensus

    def get_consensus(self,threshold,threshold_type='percent'):
        key = str(threshold)+threshold_type
        if key in self.consensuses.keys():
            return self.consensuses[key]
        return None

    def _build_stats(self):
        # Build counts for each column in the alignment                                                                               
        self.col_counts = []
        composite_count = 0
        height = len(self._alns)

        for col_num in range(len(self.composite)): # process each column                                                          
            char_counts = enumerate_column(list(self._alns.values()), col_num)
            self.col_counts.append(char_counts)
            if '-' in char_counts:
                count_base = height - char_counts['-']
            else:
                count_base = height
            composite_count += count_base
            self.col_counts[col_num]['total'] = count_base # Need this value for                                                      
        self.composite_score = float(composite_count)/height/len(self.composite)

        # Get sequence lengths from alignments                                                                                        
        self.sequence_lengths = []
        for aln in list(self._alns.values()):
            self.sequence_lengths.append(len(aln.replace('-','')))

        # Calculate num columns >= threshold vs threshold                                                                             
        self._construct_conserved_counts()

        self._reset_stats = 0

    # For each threshold value, calculates the number of columns that meet or surpass the value                                                        
    def _construct_conserved_counts(self):
        raw_count_map = {}
        # Generate map from percent conservation to number of columns with that level of conservation                                                  
        for col_num in range(len(self.composite)):
            if self.col_counts[col_num]['total'] in raw_count_map.keys():
                raw_count_map[self.col_counts[col_num]['total']] += 1
            else:
                raw_count_map[self.col_counts[col_num]['total']] = 1

        conserved_cols = []

        num_cumulative_cols = 0
        cumulative_cols = []

        summed_conservation = 0
        conservation_score = []

        #median_seq_len = median(self.sequence_lengths)
        #for conservation in range(len(self._alns)+1):
        for conservation_count in range(len(self._alns),-1,-1):
            # percent of sequences that align with a column = number of sequences that align / total number of sequences
            conservation_lvl = conservation_count/len(self._alns)

            # For this %conservation, the total number of columns that that level of conservation
            if conservation_count not in raw_count_map.keys():
                raw_count_map[conservation_count] = 0
            conserved_cols.append((conservation_lvl,raw_count_map[conservation_count]))

            # For this %conservation, the total number of columns with that level of conservation or higher
            num_cumulative_cols += raw_count_map[conservation_count]
            cumulative_cols.append((conservation_lvl,num_cumulative_cols))

            # The conservation score at threshold equal to this %conservation
            summed_conservation += conservation_lvl*raw_count_map[conservation_count]
            #conservation_score.append((conservation_lvl,summed_conservation/median_seq_len))
            conservation_score.append((conservation_lvl,summed_conservation/num_cumulative_cols))

        self._conserved_columns = list(reversed(conserved_cols))
        self._cumulative_conserved_columns = list(reversed(cumulative_cols))
        #self._score_by_threshold = list(reversed(conservation_score))
        self._score_by_threshold = list(reversed(conservation_score))

    def score_by_threshold(self):
        if self._reset_stats:
            self._build_stats()
        return self._score_by_threshold.copy()

    def percent_columns_conserved_by_threshold(self):
        if self._reset_stats:
            self._build_stats()
        percent_cols_conserved = self._cumulative_conserved_columns.copy()
        # Calculate %columns conserved from number of columns conserved/median sequence length (i.e. the benchmark number of columns)
        percent_cols_conserved = []
        for conserved_tuple in self._cumulative_conserved_columns:
            percent_cols_conserved.append((conserved_tuple[0],conserved_tuple[1]/median(self.sequence_lengths)))
        return percent_cols_conserved

    # Find k such that average conservation > k% AND columns above threshold/median seq length > k%
    def find_conservation_boundary(self):
        num_seqs = len(self._alns)
        incomplete = True
        score_by_threshold = self.score_by_threshold()
        median_seq_len = median(self.sequence_lengths)

        # Store %conserved characters from visited thresholds
        conserved_chars = {}

        # Initialize threshold to 50% converted to a count based on the number of sequences in the MSA
        thresh_count = int(0.5*num_seqs+0.5)
        above_thresh_count = num_seqs
        below_thresh_count = 0

        # TODO: Either work up from threshold of 0, or implement optimal threshold adjustment for search
        checked_counts = []
        while incomplete:
            #num_columns_above = score_by_threshold[thresh_count][1]
            num_columns_above = self._cumulative_conserved_columns[thresh_count][1]
            percent_columns = num_columns_above/median_seq_len
            sum_conservation = 0

            avg_percent_conservation = self._score_by_threshold[thresh_count][1]

            conserved_chars[thresh_count] = percent_columns

            #print("%d %.2f %.2f" % (thresh_count, avg_percent_conservation, percent_columns))

            # Determine whether the search has been exhausted
            checked_counts.append(thresh_count)
            if avg_percent_conservation == percent_columns:
                incomplete = False
                threshold = thresh_count/num_seqs
                k = avg_percent_conservation
            elif avg_percent_conservation > percent_columns:
                # Need to reduce threshold
                if thresh_count - below_thresh_count == 1:
                    # Found k
                    incomplete = False
                    threshold = thresh_count/num_seqs
                    k = avg_percent_conservation
                else:
                    above_thresh_count = thresh_count
                    thresh_count = int((thresh_count + below_thresh_count)/2)
            else: # avg_percent_conservation < percent_columns
                # Need to increase threshold
                if above_thresh_count - thresh_count == 1:
                    incomplete = False
                    threshold = above_thresh_count/num_seqs
                    k = self._score_by_threshold[above_thresh_count][1]
                else:
                    below_thresh_count = thresh_count
                    thresh_count = int((thresh_count + above_thresh_count)/2)

        return threshold,k

    def build_consensus(self, threshold_ratio, threshold_type='percent'):
        '''
        Iterates through each column, extracts respective values, counts its
        abundance and determines the consensus character given various neural
        logic clauses.
        Currently this might break if an expanded character set (of node types beyond A,C,T) is used
        '''

        height = len(self._alns) # number of entries comprising alignment
        width = len(self.composite) # all alignments are the same length

        if threshold_type == 'sqrt':
            threshold_count = threshold_ratio*math.sqrt(height) # consensus threshold
        else: # Default to percent
            threshold_count = threshold_ratio*height # consensus threshold

        consensus = None # actual generated consensus sequence
        consensus_counts = []
        #score_over_threshold = 0

        median_seq_length = median(self.sequence_lengths)

        # store consensus as list (initially) so characters can be easily
        # placed easily; akin to array indexing.
        consensus = ['-'] * width
        consensus_counts = [0] * width
        count_sum = 0 # for consensus score
        consensus_length = 0

        if self._reset_stats:
            self._build_stats()

        # Determine which character goes in each column for the consensus
        for col_num in range(width): # process each column
            char_counts = self.col_counts[col_num].copy()
            del char_counts['total'] # delete 'total' from column copy
            if '-' in char_counts:
                count_base = height - char_counts['-']
            else:
                count_base = height
            consensus_counts[col_num] = count_base

            if count_base >= threshold_count:
                consensus_length += 1
                count_sum += count_base
            ############################################################
            #### Logic for when only 1 base is exclusive to a column ###
            ############################################################
            if len(char_counts) == 1:
                chars = list(char_counts.keys()) # get the bases without counts
                char = chars[0] # get the key for that character
                consensus[col_num] = char # assign consensus character
                continue
            ############################################################
            #### Logic for when only 1 base and dashes are in column ###
            ############################################################
            if len(char_counts) == 2 and '-' in char_counts:
                # since we do not know what the known base is (either A, T, C),
                # subtract 1 from the dash count and use that to guide whether
                # the known base will become part of the consensus.
                if count_base >= threshold_count:
                    del char_counts['-'] # delete dash; remaining is character
                    char = list(char_counts.keys())[0]
                    consensus[col_num] = char # assign consensus character
                    continue
            #####################################################
            #### Logic for when 2 bases are found in a column ###
            #####################################################
            # This case should never happen, so we have an error condition
            if 'C' in char_counts and 'T' in char_counts:
                # Determine offending sequence and report in error message
                composite_is_t = self.composite[col_num] == 'T'
                offender = None
                for alignment_name in list(self._alns.keys()):
                    #print(str(alignment_ind)+' '+str(self.composite_alignments[max(alignment_ind-1,0)][0])+'\n'+str(self.composite_alignments[max(alignment_ind-1,0)][1]))
                    alignment = self._alns[alignment_name]
                    if (alignment[col_num] == 'T') ^ composite_is_t:
                        offender = alignment
                        offender_name = alignment_name
                        print('offender name: '+str(offender_name))
                if offender == None:
                    raise Exception('C and T nodes found in the same column\nPWM: '+str(self.pwm.pwm)+'\nOffender undetected')
                else:
                    raise Exception('C and T nodes found in the same column')
                    #raise Exception('C and T nodes found in the same column\nPWM: '+str(self.pwm.pwm)+'\nComposite: '+self.composite+'\nErrSeq: '+offender+'\n'+str(self.composite_alignments[max(offender_ind-1,0)][0])+'\n'+str(self.composite_alignments[max(offender_ind-1,0)][1])+'\n'+self.queries[offender_ind].seq)

                count_C, count_T = char_counts['C'], char_counts['T']
                if count_C >= count_T and count_C >= self.threshold:
                    consensus[col_num] = 'C' # select C over T
                elif count_T >= count_C and count_T >= self.threshold:
                    consensus[col_num] = 'T' # select T over C
            # choosing the better of A or C
            if 'C' in char_counts and 'A' in char_counts:
                # we get the counts for both C and A, and contrast their
                # respective scores; assigning to the consensus whichever is
                # not only the largest score but also exceed the threshold.
                if char_counts['A'] >= threshold_count:
                    consensus[col_num] = 'A' # There are enough A's
                elif char_counts['A'] + char_counts['C'] >= threshold_count:
                    consensus[col_num] = 'C' # There aren't enough A's, but there are enough characters

        consensus = ''.join(consensus)
        average_conservation = float(count_sum)/height/consensus_length

        #conservation_score = float(count_sum)/height/median_seq_length
        consensus_counts = consensus_counts
        #consensus_object = ConsensusSequence(consensus,consensus_score,consensus_counts)
        consensus_object = ConsensusSequence(consensus,threshold_ratio,threshold_type,height,average_conservation)
        return consensus_object

def enumerate_column(alns, num):
    '''
    Parses a specific column and counts the number of times a specific
    value is found in that respective column.
    '''

    char_counts = dict(Counter(''.join(aln[num] for aln in alns)))
    return char_counts
