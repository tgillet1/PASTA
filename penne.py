''' 
Performs multiple sequence alignment (MSA) given a list of
user-provided input sequences as a FASTA file.
'''

from parameter import MultipleAlignmentArgumentValidator, MultipleAlignmentCommandParser, InputWrapperState
from pairwise import PairwiseDriver
from msa import MultipleSequenceDriver, ConsensusFilterFactory
from random import shuffle
os.system("taskset -p 0xff %d" % os.getpid())

version = 1.0

def run_msa(queries, input_state):
    args = input_state.args
    if args['subsample'] != 1:
        if args['subsample'] < 1:
            subsample_size = int(args['subsample']*len(queries))
        elif args['subsample'] > 1:
            subsample_size = int(min(len(queries),args['subsample']))
            
        if args['random_subset']:
            # Randomize order and take first _subsample_size_ sequences, then reorder
            query_ids = [x for x in range(len(queries))]
            shuffle(query_ids)
            query_ids = sorted(query_ids[0:subsample_size])
            queries = list([queries[i] for i in query_ids])
        else:
            if args['subsample_start'] < 1:
                start = int(args['subsample_start']*len(queries))
            else:
                start = int(args['subsample_start'])
            queries = list([queries[i] for i in range(start,start+subsample_size)])

    if args['random_order']:
        shuffle(queries)

    # Generate MSA
    msa_driver = MultipleSequenceDriver(queries, input_state)
    msa = msa_driver.get_msa()

    # Build resultant consensus
    consensus_object = msa.build_consensus(args['threshold'],args['type'])

    # Write MSA and consensus to file
    #consensus_fact = ConsensusFilterFactory(msa_driver,consensus_object)
    consensus_fact = ConsensusFilterFactory(msa,consensus_object)
    consensus_fact.write(fname=args['build'])

    #consensus_fact = ConsensusFilterFactory(driver.alns, driver.composite, args['threshold'], args['type'])
    #consensus_fact.build_consensus()
    #consensus_fact.write(fname=args['build']) # write analysis
    
if __name__ == '__main__':
    try:
        args = MultipleAlignmentCommandParser().parse_args()
        MultipleAlignmentArgumentValidator(args) # test all arguments are correct

        print('penne - v.' + str(version) + '\n=================')
        input_state = InputWrapperState(args)
        input_state.assign_matrix() # parse in-built or custom matrix
        targets = input_state.parse_fasta(input_state.fname) # next, parse fasta file
        
        run_msa(queries, input_state)
            
    except (IOError, KeyboardInterrupt, IndexError) as e:
        print(str(e)+'\n')
