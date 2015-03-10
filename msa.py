from pairwise import NeedlemanWunsch, PositionWeightedMatcher
import math
from sequence import NeuriteSequence, MultipleSequenceAlignment
import concurrent.futures
import sequence
from collections import Counter
# from random import shuffle
from tree import TreeLogicFactory
from buffer import XMLBuildWriter

class MultipleSequenceDriver():
    ''' 
    A high-level class to perform multiple sequence pairwise.
    If iterate is an integer the align phase will run multiple times with a pwm based on the
    prior run.
    '''
    def __init__(self, queries, input_state):
        self.queries = queries
        self.costs = input_state.get_penalties() # set costs to core
        self.submat = input_state.get_submatrix() # set submatrix to core
        self.iterate = input_state.get_args()['iterate']
        self.MSA = None # MSA must be built
        self.alns = {}
        self.composite = None # initially, no composite exists
        self.composite_score = 0
        self.alignment_file = input_state.alignment_file
        self.pwm = None
        self.node_types = input_state.node_types
        self.consensus_check_percent = .4
        self.num_workers = input_state.get_args()['num_processes']

        self.composite_alignments = []

    def get_msa(self):
        if self.MSA is None:
            self.build_composite()
            self.align()
        return self.MSA

    def build_composite(self):
        ''' 
        Takes the first 2 sequences, aligns them and renders the pairwise
        alignment a composite (formerly "pre-consensus"). Every sequence thereof is then 
        subsequently pairwise-aligned to this composite.
        '''

        print("Building composite")

        queries = self.queries
        # get the first two input sequences
        s0 = queries[0]
        s1 = queries[1]
        # pass them both into the tree--based Needleman--Wunsch algorithm.
        nw = NeedlemanWunsch(s1=s0, s2=s1, costs=self.costs, submat=self.submat, 
                                       node_types=self.node_types)
        first_align, second_align = nw.prettify()[1]
        self.composite_alignments.append([nw.align1,nw.align2])
        
        # feed respective alignments into an analysis class and get consensus.
        composite = NeuriteSequence('Composite',
                                    TreeLogicFactory(str1=first_align, 
                                                     str2=second_align).get_alignment())

#        print("Seq1: "+first_align)
#        print("Seq2: "+second_align)

#        print("Comp: "+composite)

        # since the first two sequences have been aligned, focus on all others.
        for i in range(2, len(queries)):
            curr_seq = queries[i]
            nw = NeedlemanWunsch(s1=composite, s2=curr_seq, 
                                 costs=self.costs, submat=self.submat, 
                                 node_types=self.node_types)

            align_sA, align_sB = nw.get_alignment()

            self.composite_alignments.append([align_sA,align_sB])
            #print(align_sA)
            #print(align_sB)

            composite = NeuriteSequence('Composite',
                                        TreeLogicFactory(str1=align_sA, 
                                                         str2=align_sB).get_alignment())
            self.composite = composite # of type NeuriteSequence        

    def align(self):
        ''' 
        A composite alignment is a manifestation of pairwise alignments
        between sequence pairs. In other words, sequence (i) and (i+1) are
        aligned and the resultant alignment (composite) is saved. Next,
        sequence (i+2) is aligned against this composite and the resultant
        alignment updates the composite. This logic is repeated for all
        queries. This function re-aligns the input queries back onto the
        composite so that the actual consensus sequence can be derived.
        If the 'iterate' parameter is set, the last part will repeat using
        a pwm produced using the previous run.
        '''
        if not self.composite:
            raise IOError('Composite alignment required.')
        if self.alignment_file:
            align_handle = open(self.alignment_file, 'w')
            
        name_lengths = []
        iterate = self.iterate
        if iterate < 1:
            iterate_count = 1 # Iterate until %change is within the iterate threshold given in iterate

        # Initialize pwm with equal weights per position
        pwm = PositionWeightedMatrix(self.composite,list([{'total':1,char:1} for char in self.composite.seq]))
        continue_iter = True
        iter_count = 0
        change = 1
        consensus_conservation_score = 0
        consensus_length = 0
        #total_score_over_threshold = 0
        prev_consensus_score = 0
        prev_consensus_length = 0
        # Iternate until change in consensus strength is sufficiently small or some number of iterations have been run
        while (iterate < 1 and (change > iterate or consensus_length > prev_consensus_length)) or iterate >= 1 and iter_count < iterate:
#        while (iterate < 1 and change > iterate) or (iterate >= 1 and iter_count < iterate):
            iter_count += 1

#        for curr_it in range(iterate_count):
            self.alns = {} # New alignments at each iteration, final iteration gives final alignment
            executor = concurrent.futures.ProcessPoolExecutor(self.num_workers)
            try:
                # Align each query with composite
                for curr_seq in self.queries:
                    # Get name lengths for position aligning names in alignment file output
                    if len(name_lengths) == 0:
                        name_lengths.append(len(curr_seq.name))

                    # Setting 'consensus=2' tells NW that s2 is the consensus and will prevent gaps from appearing in s1 alignemtn

                    #pw_matcher = PositionWeightedMatcher(sequence=curr_seq, pwm=pwm, 
                    #                    costs=self.costs, node_types=self.node_types)
                    f = executor.submit(PositionWeightedMatcher,sequence=curr_seq, pwm=pwm, 
                                        costs=self.costs, node_types=self.node_types)
                    f.add_done_callback(self._msa_callback)
                executor.shutdown()
            except KeyboardInterrupt:
                executor.shutdown()
                print("MSA Iteration interupted, stopping program")
                change = iterate + 1
                iter_count = iterate

            # Get the PWM given the alignments for the next iteration
            pwm = self.get_pwm()

            # Calculate change between previous and current alignment (even if not required)
            prev_consensus_score = consensus_conservation_score
            prev_consensus_length = consensus_length

            self.MSA = MultipleSequenceAlignment(self.composite.seq,self.alns,node_types=self.node_types)
            consensus_check_obj = self.build_consensus(threshold_ratio=self.consensus_check_percent)
            consensus_conservation_score = consensus_check_obj.average_conservation
            consensus_length = consensus_check_obj.get_consensus_length()
            change = consensus_conservation_score - prev_consensus_score

            print("MSA Iteration "+str(iter_count)+"; consensus length: "+str(consensus_length)+"; consensus score: "+str(consensus_conservation_score))

        # Write alignments to file
        if self.alignment_file:
            total_space = max(12,max(name_lengths))+1
            align_handle.write(('Composite'+' '*(total_space-12))+self.composite.seq+'\n') # write header
            for curr_seq in self.queries:
                align_handle.write(curr_seq.name+(' '*(total_space-len(curr_seq.name)))+(''.join(self.alns[curr_seq.name]))+'\n') # write header
            align_handle.close()

        # Create MSA object from composite and final alignment
        self.MSA = MultipleSequenceAlignment(self.composite.seq,self.alns,self.node_types)

    def _msa_callback(self, return_val):
        cbexcept = return_val.exception()
        if cbexcept is not None:
            print("Error: "+str(cbexcept))
        pw_matcher = return_val.result()
        align_sA = pw_matcher.get_alignment()[0]
        #self.alns.append(align_sA)
        self.alns[pw_matcher.seq1.name] = align_sA


    # Generate a position weighted matrix of the form: array of {character:weight}
    def get_pwm(self):        
        # Determine composite score now that alignments exist, and generate pwm
        col_counts = []
        composite_count = 0
        height = len(self.queries)

        for col_num in range(len(self.composite.seq)): # process each column
            char_counts = enumerate_column(list(self.alns.values()), col_num)
            col_counts.append(char_counts)
            if '-' in char_counts:
                count_base = height - char_counts['-']
            else:
                count_base = height
            composite_count += count_base
        self.composite_score = float(composite_count)/height/len(self.composite.seq)

        composite_node_types = NeuriteSequence('PWM Node Types',''.join(list([self.node_types[char] for char in self.composite.seq])))
        self.pwm = PositionWeightedMatrix(composite_node_types,col_counts)

        return self.pwm

    def build_consensus(self, threshold_ratio, threshold_type='percent'):
        return self.MSA.build_consensus(threshold_ratio,threshold_type)

class PositionWeightedMatrix():
    '''
    Contains a string of node types (A,C,T) along with an array of dictionaries giving the weight of a position and
    an associated character.
    '''
    def __init__(self, node_type_sequence, pwm):
        self.node_type_sequence = node_type_sequence
        self.pwm = pwm
        # Calculate total weight of each column
        for pos in range(len(self.pwm)):
            self.pwm[pos]['total'] = sum(self.pwm[pos][key] for key in self.pwm[pos] if not key=='total' and not key=='-')

def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]

class ConsensusFilterFactory():
    ''' 
    When multiple sequence alignment is complete, alignments per query against
    the composite are produced. Each alignment must then be parsed given
    neural logic rules, enabling derivation of a sound consensus sequence.
    '''
    
#    def __init__(self, alignments, composite, threshold_ratio, threshold_type):
#        self.alns = alignments # matrix representing alignments
#        self.threshold = threshold_ratio
    def __init__(self, msa, consensus_object):
        self.msa = msa
        self.consensus_obj = consensus_object
            
    def write(self, fname):
        ''' 
        Saves MSA analysis as an XML file.
        @param fname: Output filename
        '''
        #w = XMLBuildWriter(fname, self.msa_obj, self.consensus_obj)
        w = XMLBuildWriter(fname, self.msa, self.consensus_obj)
        w.write() # write the XML data-structure to the disk
    

def enumerate_column(alns, num):
    '''
    Parses a specific column and counts the number of times a specific
    value is found in that respective column.
    '''

    char_counts = dict(Counter(''.join(aln[num] for aln in alns)))
    return char_counts

