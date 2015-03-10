''' 
Determines the amount of time taken for aligning a set of sequences of a given size
'''

from parameter import AlignmentArgumentValidator, AlignmentCommandParser, InputWrapperState
from pairwise import PairwiseDriver
from timeit import Timer
import time
import csv

version = 0.1
    
if __name__ == '__main__':
    try:
        #args = AlignmentCommandParser().parse_args()
        #AlignmentArgumentValidator(args) # test all arguments are correct

        gen_args = {'f': 'A', 'f2': None, 'a': None, 'gap': -1, 'custom': '/home/tgillett/PairwiseAlignment/identity.tab', 'n': 1}

        #print('spaghetti - v.' + str(version) + '\n=================')

        matrix = [['Length','G0','G1','G3']]
        fasta_dir = '/home/tgillett/fastaFromSwc/random/Mapping2/Method3/'

        input_state = InputWrapperState(gen_args)
        input_state.assign_matrix() # parse in-built or custom matrix

        for seqlen in [16,18,20,22,24,26,29,32,35,39,43,47,52,57,63,69,76,84,92,101,111,122,134,147,162,178,196,216,238,262,288,317,349,384,400,422,432,464,504,544,588,635]:
            queries = input_state.parse_fasta(fasta_dir+'RandASTL-Length'+str(seqlen)+'_1000Seqs.fasta')
            targets = input_state.parse_fasta(fasta_dir+'RandBSTL-Length'+str(seqlen)+'_1000Seqs.fasta')
            times = [seqlen]
            for gap_open in [0,-1,-3]:
                args = gen_args

                args['gapopen'] = gap_open
                args['o'] = '/home/tgillett/PairwiseAlignment/Baseline/GapOpen'+str(gap_open)+'_1000/BaselineScoreSTL-Lengths'+str(seqlen)+'_'+str(seqlen)+'.tab'
                input_state = InputWrapperState(args)

                driver = PairwiseDriver(targets, queries, input_state)
                startTime = int(round(time.time() * 1000))
                driver.start() # start only pairwise alignment
                endTime = int(round(time.time() * 1000))
                times.append(endTime-startTime)
            matrix.append(times)

        with open("AlignmentTimes.csv", "wb") as f:
            writer = csv.writer(f)
            writer.writerows(matrix)

    except (IOError, KeyboardInterrupt, IndexError) as e:
        print(str(e)+'\n')
