''' 
Performs both local sequence alignment given a list of
user-provided input sequences as a FASTA file.
'''

from parameter import PairwiseAlignmentArgumentValidator, PairwiseAlignmentCommandParser, InputWrapperState
from pairwise import PairwiseDriver
from msa import MultipleSequenceDriver, ConsensusFilterFactory
from random import shuffle
from score_converter import ScoreConversionDriver
import os
os.system("taskset -p 0xff %d" % os.getpid())

version = 1.0

if __name__ == '__main__':
    try:
        args = PairwiseAlignmentCommandParser().parse_args()
        PairwiseAlignmentArgumentValidator(args) # test all arguments are correct

        print('spaghetti - v.' + str(version) + '\n=================')
        input_state = InputWrapperState(args)
        input_state.assign_matrix() # parse in-built or custom matrix
        queries = input_state.parse_fasta(input_state.fname) # next, parse fasta file
        if input_state.fname2 is None:
            targets = queries
        else:
            targets = input_state.parse_fasta(input_state.fname2)

        if input_state.get_args()['mode'] == 'align':
            driver = PairwiseDriver(targets, queries, input_state)
            driver.start() # start only pairwise alignment
        elif input_state.get_args()['mode'] == 'calculate_percharacter':
            if input_state.get_args()['rawscore_file'] is None:
                raise IOError('While using mode "calculate_percharacter", a rawscore_file must be provided.')
            driver = ScoreConversionDriver(input_state,queries,targets)
            driver.start()
            
    except (IOError, KeyboardInterrupt, IndexError) as e:
        print(str(e)+'\n')
