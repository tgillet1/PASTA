''' 
Runs some finer statistical analysis of MSAs. The initial version is meant to give a value
for the weighted amount of sequence material between characters in the consensus to represent
the variability of those regions. Consensus characters with no gaps between represent a 
level of variability no greater than that defined by the consensus threshold.
'''

from parameter import ConsensusStatsCommandParser
from buffer import XMLBuildReader, status_message
import pairwise
from msa import MultipleSequenceDriver, ConsensusFilterFactory
from sequence import NeuriteSequence

version = 0.1

def calculate_conservation(msa,consensus):
    raw_consensus = consensus.raw_consensus
    alignments = msa.get_alignments()

    contributions = {}
    for name in list(alignments.keys()):
        alignment = alignments[name]
        contributions[name] = 1/len(alignment.replace('-',''))

    current_var_sum = 0
    variability_sums = []
    conservation = []

    ## Calculate conservation at and variability before each character

    # Go through each position in the alignment (based on the composite or gapped consensus)
    for position in range(len(msa.composite)):
        # If the current position is a gap, add contributions to the current variability sum
        if raw_consensus[position] == '-':
            # Go through each alignment to check whether it has a character at the current position
            for name in list(alignments.keys()):
                alignment = alignments[name]
                # If it does have a character (not a gap '-'), add its variability contribution
                if not alignment[position] == '-':
                    current_var_sum += contributions[name]
        # The current position is a consensus character, so store the variability sum preceeding it
        else:
            # if the position within the ungapped consensus is needed, just use the current length of variability_sums
            variability_sums.append(current_var_sum)
            current_var_sum = 0

            # determine the current position's conservation level
            conservation_sum = 0
            for alignment in list(alignments.values()):
                # If it does have a character (not a gap '-'), add its conservation contribution
                if not alignment[position] == '-':
                    conservation_sum += 1
            conservation.append(conservation_sum/len(alignments))

    return conservation,variability_sums

#TODO - A stack is going to be required for this. Look at Java code to recreate (shouldn't be too hard)
def generate_consensus_newick(msa,consensus,filename):
    newick_string = ""
    cons_str = consensus.consensus
    conservation,variability = calculate_conservation(msa,consensus)
    a_stack = [(1,0)]
    print(cons_str)
    for position in range(len(cons_str)):
        #print(str(position)+' '+cons_str[position])
        #print(str(a_stack))
        if cons_str[position] == 'A':
            a_stack.append((0,0))
            node_string = ')'
        elif cons_str[position] == 'C':
            node_string = ',:1)'
            a_val = a_stack.pop()
            a_stack.append((a_val[0],a_val[1]+1))
        elif cons_str[position] == 'T':
            node_string = '(:1,:1)'
            not_done = 1
            while not_done and len(a_stack) > 0:
                a_val = a_stack.pop()
                if a_val[0] == 0:
                    node_string = ','+a_val[1]*'('+node_string
                    a_stack.append((1,1))
                    not_done = 0
                else:
                    node_string = a_val[1]*'('+node_string
            
        newick_string = node_string + cons_str[position] + '-' + str(round(conservation[position],3)) + ':' + str(round(variability[position],3)) + newick_string

    newick_string += ';'

    status_message('Generating newick string','OK')
    handle = open(args['newick'],'w')
    handle.write(newick_string)
    handle.close()
    return newick_string

def generate_score_and_conserved_chars_file(msa,filename):
    score_by_threshold = msa.score_by_threshold()
    percent_cols_conserved = msa.percent_columns_conserved_by_threshold()
    handle = open(args['scores'],'w')
    for i in range(len(score_by_threshold)):
        handle.write("%.3f,%.3f,%.3f\n" % (score_by_threshold[i][0],score_by_threshold[i][1],percent_cols_conserved[i][1]))
        #print("%.3f  %.3f  %.3f" % (score_by_threshold[i][0],score_by_threshold[i][1],percent_cols_conserved[i][1]))
    handle.close()
    print("Wrote scores by threshold to "+filename)

if __name__ == '__main__':
    try:
        args = ConsensusStatsCommandParser().parse_args()
#        ConsensusStatsArgumentValidator(args) # test all arguments are correct

        print('ditalini - v.' + str(version) + '\n=============')

        msa = XMLBuildReader(args['build']).parse()

        if args['threshold']:
            threshold = args['threshold']
            msa.add_consensus(threshold,msa.build_consensus(float(threshold)))
            print('Consensus at threshold ' + threshold + ': ' + msa.get_consensus(threshold).consensus)

        if args['newick']:
            if args['threshold']:
                generate_consensus_newick(msa,msa.get_consensus(float(args['threshold'])),args['newick'])
            else:
                generate_consensus_newick(msa,list(msa.consensuses.values())[0],args['newick'])

        if args['scores']:
            generate_score_and_conserved_chars_file(msa,args['scores'])                

        if args['k']:
            threshold,k = msa.find_conservation_boundary()
            print("At threshold %.2f, average conservation and proportion of conserved characters are both %.2f%%" % (threshold, k*100))

        if args['newbuild']:
#            consensus_object = msa.build_consensus(args['threshold'],args['type'])

            # Write MSA and consensus to file
            consensus_fact = ConsensusFilterFactory(msa,msa.get_consensus(threshold))
            consensus_fact.write(fname=args['newbuild'])

        status_message('Consensus statistics computation complete ', 'OK')

    except (IOError, KeyboardInterrupt, IndexError) as e:
        print(str(e))
