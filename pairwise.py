from matrix import StateMatrix, DirectionalMatrixWrapper
from buffer import status_message
import concurrent.futures
import traceback
import sys
from sequence import NeuriteSequence
from itertools import repeat
import os
# from random import shuffle
os.system("taskset -p 0xff %d" % os.getpid())

def create_ta_dictionary(seq,node_types,submatrix=None,gap_cost=-1):
	'''
	Creates a dictionary for a given tree sequence linking each T node to 
	an associated A node.
	'''
	if submatrix is None:
		submatrix = {}

	for nodetype in node_types:
		if not (nodetype,'-') in submatrix.keys():
			submatrix[(nodetype,'-')] = gap_cost

	# Dictionary with (key,value) pairs of (T-node index, A-node index) as well as (str(T-node index), total gap cost)
	taDict = {}
	AStack = []
	AStack.append(-1) # Begin the AStack with a sentinal value
	CostStack = []
	CostStack.append(0)
	# Visit each character in the sequence
	for index in range(0,len(seq)):
		# Add the cost of the current node to the current gap cost register
		CostStack[-1] += submatrix[seq[index],'-']

		# If the current node is an A-node
		if seq[index] in node_types['A']:
			# push it onto the top of the stack
			AStack.append(index)
			# push on a new gap cost register
			CostStack.append(0)
			# The cost of an A-node at the front of a gap is added on in
			# NeedlemanWunsch.calculate_gap when determining whether the A is gapped or is
			# matched to a C-node

		# If the current node is a T-node
		elif seq[index] in node_types['T']:
			# Set the dictionary to contain the (T-node index, associated A-node pair)
			taDict[index] = AStack.pop()
			# Get the cost of the gap between the T and A nodes
			currCost = CostStack.pop()
			# Add that cost to the dictionary
			taDict[str(index)] = currCost
			# Add the cost to the enclosing AT pair
			if (len(CostStack) > 0):
				CostStack[-1] += currCost
	return taDict


class PairwiseDriver():
    '''
    Executes the local alignment application
    '''
    def __init__(self, targets, queries, input_state=None, store_pairwise=False, score_type='alignment_score', **kwargs):
        self.targets = targets # core operates given targets and queries
        self.queries = queries
        self.query_map = {query.name:query for query in queries}
        self.num_complete = 0 # len(self.priorCompletions) # for how many sequences have been aligned
        self.store_pairwise = store_pairwise
        if store_pairwise:
            self.score_dict = {}
        self.score_type = score_type # Can also be 'num_gaps'
        acceptable_score_types = ['alignment_score','num_gaps','per_character_score']
        if not score_type in acceptable_score_types:
            raise IOError('Unknown score type "'+score_type+'". Acceptable values: '+ acceptable_score_types )

        if input_state is None:
            input_state = InputStateWrapper(kwargs)

        self.costs = input_state.get_penalties() # set costs to core
        self.submat = input_state.get_submatrix() # set submatrix to core
        scoremat_outfile = input_state.get_args()['output']
        
        # Set openMode to append if some targets have already been run and completed
        if self.num_complete > 0:
            openMode = 'a'
        else:
            openMode = 'w'
            
        self.scorehandle = None
        if not scoremat_outfile is None:
            self.scorehandle = open(scoremat_outfile, openMode) # output file

        self.alignhandle = None
        if input_state.get_args()['alignment_file'] is not None:
            self.alignhandle = open(input_state.get_args()['alignment_file'], openMode) # alignments file
            
        self.num_workers = input_state.get_args()['num_processes']
        # Get node type lists
        self.node_types = input_state.get_node_types()

        self.debug = 0

    def set_debug(val):
        self.debug = val

    # Initialize the core given query sequences and input arguments
    def start(self):
        status_message('Pairwise alignment running', 'please wait')
        executor = concurrent.futures.ProcessPoolExecutor(self.num_workers)
        try:
            with concurrent.futures.ProcessPoolExecutor(self.num_workers) as executor:
                for result in executor.map(_aligner, self.targets, repeat(self.queries), repeat(self.costs), 
                                           repeat(self.submat), repeat(self.node_types)):
                    self._process_result(result)
#            for target in self.targets: # per fasta, create a concurrent job, f.
#                print(target.name)
#                f = executor.submit(_aligner, target, 
#                                    self.queries, self.costs, 
#                                    self.submat, self.node_types)
#                f.add_done_callback(self._callback)
            executor.shutdown()
            self.close_output_buffers()
            status_message('Analysis complete', 'OK')
        except KeyboardInterrupt:
            executor.shutdown()

    # Close all I/O buffers such as file handles
    def close_output_buffers(self):
        if self.alignhandle is not None:
            self.alignhandle.close()
        if self.scorehandle is not None:
            self.scorehandle.close()

    # Get the headers, i.e. top-most row for the score matrix
    def _create_header(self, results):
        h = '\t' +'\t'.join([h[-1] for h in results])
        self.scorehandle.write(h + '\n')
        self.scorehandle.flush()

    def _process_result(self, result):
        target, results = result
        #results = sorted(results, key=lambda x: x[-1]) # sort by query (last item)

        # determine scores
        if self.score_type == 'alignment_score':
            scores = [result[0] for result in results]
        elif self.score_type == 'num_gaps':
            scores = [result[1][0].count('-')+result[1][1].count('-') for result in results]
        elif self.score_type == 'per_character_score':
            scores = [(result[0] + abs(self.costs['gapopen']) + 
       	                  abs(self.costs['gap'] * (len(target.seq)-len(self.query_map[result[2]].seq)))) / 
                        min(len(target.seq),len(self.query_map[result[2]].seq)) 
                      for result in results]

	# save scores
        if self.store_pairwise:
            self.score_dict[target] = scores

        # write scores to output file
        if self.scorehandle is not None:
            if self.num_complete == 0: # for the first result, write headers
                self._create_header(results)
            scores = '\t'.join([str(s) for s in scores])
            self.scorehandle.write(target + '\t' + scores + '\n')
            self.scorehandle.flush()

        # also save actual alignment string
        if self.alignhandle is not None:
            for r in results:
                if r[1] is not None:
                    align_target, align_query = r[1]
                    out_str = target + '\t' + r[-1] +'\t'+ align_target +'\t'+ align_query
                    self.alignhandle.write(out_str + '\n')
                    self.alignhandle.flush()
        self.num_complete += 1
        if self.debug >= 1:
            print(' --> ' + target + ' [OK] '+str(self.num_complete) +
                ' of '+ str(len(self.targets))) # print progress
        
        
    # Callback function once a thread is complete
    def _callback(self, return_val):
        cbexcept = return_val.exception()
        if cbexcept is not None:
            print("Err1: "+str(cbexcept))
            #print("Targets: "+str(list([x.seq for x in self.targets])))
            #print("Queries: "+str(list([x.seq for x in self.queries])))
        res = return_val.result() # get result once thread is complete

        self._process_result(res)

#        target, results = res
#        #results = sorted(results, key=lambda x: x[-1]) # sort by query (last item)
#
#        # determine scores
#        if self.score_type == 'alignment_score':
#            scores = [result[0] for result in results]
#        elif self.score_type == 'num_gaps':
#            scores = [result[1][0].count('-')+result[1][1].count('-') for result in results]
#        elif self.score_type == 'per_character_score':
#            scores = [(result[0] + abs(self.costs['gapopen']) + 
#       	                  abs(self.costs['gap'] * (len(target.seq)-len(self.query_map[result[2]].seq)))) / 
#                        min(len(target.seq),len(self.query_map[result[2]].seq)) 
#                      for result in results]
#
#	# save scores
#        if self.store_pairwise:
#            self.score_dict[target] = scores
#
#        # write scores to output file
#        if self.scorehandle is not None:
#            if self.num_complete == 0: # for the first result, write headers
#                self._create_header(results)
#            scores = '\t'.join([str(s) for s in scores])
#            self.scorehandle.write(target + '\t' + scores + '\n')
#            self.scorehandle.flush()
#
#        # also save actual alignment string
#        if self.alignhandle is not None:
#            for r in results:
#                if r[1] is not None:
#                    align_target, align_query = r[1]
#                    out_str = target + '\t' + r[-1] +'\t'+ align_target +'\t'+ align_query
#                    self.alignhandle.write(out_str + '\n')
#                    self.alignhandle.flush()
#        self.num_complete += 1
#        if self.debug >= 1:
#            print(' --> ' + target + ' [OK] '+str(self.num_complete) +
#                ' of '+ str(len(self.targets))) # print progress

    def get_score_matrix(self):
        score_mat = []
        for target in self.targets:
            if target.name in self.score_dict:
                score_mat.append(self.score_dict[target.name])
        return score_mat
        
# Maps each query sequence against a set of targets (itself)
def _aligner(target, queries, costs, submat, node_types):
    results = [] # K => target, V => aligned queries 
    # get the gap and substitution matrix
    for query in queries:
        NW = NeedlemanWunsch(target, query, costs, submat, node_types)
        output = NW.prettify()
        results.append(output)
    return target.name, results

# Implementation of global alignment - based on Needleman-Wunsch
class NeedlemanWunsch():
	def __init__(self, s1, s2, costs, submat, node_types, composite=0):
		self.seq1 = s1 # sequence 1
		self.seq2 = s2 # sequence 2
		self.costs = costs # dictionary of all costs (i.e. penalties)
		self.submat = submat # substitution matrix
		self.create_node_types(node_types)
		self.create_residue_specific_gapcost()
		self.TADict1 = create_ta_dictionary(s1.seq, node_types, submat, costs['gap'])
		self.TADict2 = create_ta_dictionary(s2.seq, node_types, submat, costs['gap'])
		self.scoreMat = None # references score matrix
		self.directionMat = None # references diag(0),left(1),up(2) matrix
		self.leftMat = None # references diag(0),left(1),up(2) matrix
		self.upMat = None # references diag(0),left(1),up(2) matrix
		self.backPos = {} # the backtrace position from one position to its prior (contains integer pairs)
		self.align1 = '' # alignment string for sequence 1
		self.align2 = '' # alignment string for sequence 2
		self.composite = composite
		self._aligner()
	
	def create_node_types(self,node_types):
		self.node_types = {}
		for node_type in node_types.keys():
			for residue in node_types[node_type]:
				self.node_types[residue] = node_type
	
	# Fills in all gap-residue pairs either with flipped order entry if it exists, else with the default gap cost
	# Also fills in flipped residue-residue scores
	def create_residue_specific_gapcost(self):
		newToSubmat = {}
		for pair in self.submat.keys():
			if pair[0] == '-' and pair[1] == '-':
				# do nothing, this is useless and shouldn't happen
				pass
			elif pair[0] == '-' or pair[1] == '-':
				# Note that the given residue has an associated gap cost
				if pair[0] == '-' and (pair[1],'-') not in self.submat.keys():
					newToSubmat[pair[1],'-'] = self.submat[pair]
				elif pair[1] == '-' and (pair[1],'-') not in self.submat.keys():
					newToSubmat[pair[1],'-'] = self.submat[pair]
			else:
				# Fill the residueDict so none are missed
				if (pair[0],'-') not in self.submat.keys() and ('-',pair[0]) not in self.submat.keys():
					newToSubmat[pair[0],'-'] = self.costs['gap']
					newToSubmat['-',pair[0]] = self.costs['gap']
				if (pair[1],'-') not in self.submat.keys() and ('-',pair[1]) not in self.submat.keys():
					newToSubmat[pair[1],'-'] = self.costs['gap']
					newToSubmat['-',pair[1]] = self.costs['gap']
				if (pair[1],pair[0]) not in self.submat.keys():
					newToSubmat[pair[1],pair[0]] = self.submat[pair]
		for pair in newToSubmat.keys():
			self.submat[pair] = newToSubmat[pair]
	
	# return top (highest) alignment score given sequence 1 and 2
	def get_top_score(self):
		return self.scoreMat.get_data(-1,-1)
	
	# alignment string resultant from sequence 1 and 2
	def get_alignment(self):
		return (self.align1, self.align2)
	
	# Create parser-friendly output given a NW alignment 
	def prettify(self):
		return [ self.get_top_score(), self.get_alignment(), self.seq2.name ]

	def determine_open_extend(self,i,j,m,directionM,dirScoreM,currentGapCost,gap_direction):
		if m.get_data(i,j) is None:
			return None,False
		gapScore = m.get_data(i,j) + currentGapCost
		# Add on the gap open cost if the prior position is not gapped in the same direction (gapping the same sequence)
		if dirScoreM.score.get_data(i,j) is None:#[0]): # Previous position can't gap
			gapScore += self.costs['gapopen']
			scoreExtendPair = gapScore,False
		elif directionM.get_data(i,j) is not gap_direction: # Previous position didn't choose gap
			extendGapScore = dirScoreM.score.get_data(i,j) + currentGapCost
			openGapScore = gapScore + self.costs['gapopen']
			if openGapScore >= extendGapScore: # new gap is best, go with previous position's choice
				scoreExtendPair = openGapScore,False
			else: 
				# continuation is best, override previous position choice; if current position is on path, previous position will gap
				scoreExtendPair = extendGapScore,True
		else: # Previous position did choose gap, so this is automatically a continuation
			scoreExtendPair = gapScore,True
		return scoreExtendPair
		
	# Calculates the gap cost in a given direction from a given position, which depends on the node type
	def calculate_gap(self,i,j,seq1,seq2,m,directionM,dirScoreM,TADict,gap_direction):
		isExtend = False
		a_c_match = False
		# CType: gap one
		if self.node_types[seq1[i-1]] == 'C':
			# The prior position assuming a gap (index based on m)
			gapPosi = i - 1
			gapPosj = j
			# Determine the appropriate gap score depending on whether this opens or extends a gap
			gapScore,isExtend = self.determine_open_extend(gapPosi,
														gapPosj,
														m,
														directionM,
														dirScoreM,
														get_gapcost(seq1[i-1],self.submat),
														gap_direction)

		# TType: gap until paired A
		elif self.node_types[seq1[i-1]] == 'T': 			
			# Special Case: Final position of sequence
			if i == len(seq1):
				# Check all possible previous T positions for the case in which the other sequence (seq2) is a subtree of this sequence (seq1)
				gapScore = None
				gapPosi = None
				gapPosj = None
				isExtend = False
				for k in range(1,len(seq1)-1):
					# Previous position must be a T and must match to be a candidate
					if self.node_types[seq1[k-1]] == 'T' and directionM.get_data(k,j) == 0:
#						print(str(k)+' '+str(j)+' '+str(m.score.get_data(k,j)))
#						print('\n'.join([''.join(["{:5s}".format(str(cell)) for cell in row]) for row in dirScoreM.score.data]))
						cur_score = m.get_data(k,j) + self.costs['gapopen'] - (len(seq1)-k)
						if gapScore is None or cur_score > gapScore:
							gapScore = cur_score
							gapPosi = k
				if gapScore is not None:
					gapPosj = j
				
			else: # Normal case (gap until paired A)
				# The prior position assuming a gap (index based on m)
				gapPosi = TADict[i-1]
				# Case where this is the last T; handle sentinal and get cost of front-gap
				if TADict[i-1] is -1:
					gapPosi = 0
					gapCostStart = 0
				else:
					gapCostStart = get_gapcost(seq1[gapPosi],self.submat) # Cost of the A-node that starts the gap

				gapCostMajor = TADict[str(i-1)] # Cost of the gap from the T-node up to the A-node

				# Calculate the total gap cost assuming the associated A is also gapped
				gapScoreGapFinish,isExtend = self.determine_open_extend(gapPosi,
																j,
																m,
																directionM,
																dirScoreM,
																gapCostMajor+gapCostStart,
																gap_direction)

				# If seq2 character is C-type, determine whether to match T-paired A and C, or to just gap the A
				# This will not happen if this is the last T (TADict[i-1] is -1), as the whole sequence must be gapped
				if self.node_types[seq2[j-1]] == 'C' and TADict[i-1] is not -1:
					# Calculate the total gap cost assuming the associated A-node matches a C-node
					ACScore = get_score(seq1[gapPosi],seq2[j-1],self.submat)
					if m.get_data(gapPosi,j-1) is None:
						gapScoreACFinish = None
					else:
						gapScoreACFinish = m.get_data(gapPosi,j-1) + gapCostMajor + ACScore + self.costs['gapopen']

					# Determine which gap produces a higher overall score, and use that for this position's gap score
					if gapScoreACFinish is not None and (gapScoreGapFinish is None or gapScoreACFinish >= gapScoreGapFinish):
						# Best score achieved if gap finishes with A-C match
						gapScore = gapScoreACFinish
						gapPosj = j-1
						isExtend = False
						a_c_match = True
					elif gapScoreGapFinish is not None:
						# Best score achieved if gap includes the associated A
						gapScore = gapScoreGapFinish
						gapPosj = j
					else:
						# Gap cannot be completed
						gapScore = None
						gapPosj = j
						isExtend = False
				# If seq2 character is not a C, then use the total gap score assuming the A is also gapped
				else:
					gapScore = gapScoreGapFinish
					gapPosj = j
		else: # AType, no gapping allowed
			gapScore = None # original was set to None
			gapPosi = None
			gapPosj = None 
		dirScoreM.score.set_data(i,j, gapScore)
		dirScoreM.extend_flag.set_data(i,j, isExtend)
		dirScoreM.a_c_match.set_data(i,j, a_c_match)
		return [gapScore, gapPosi, gapPosj]

	# Execute alignment
	def _aligner(self):
		l1, l2 = len(self.seq1.seq), len(self.seq2.seq)	
		self.scoreMat = StateMatrix(l1+1, l2+1) # create matrix for storing counts
		self.directionMat = StateMatrix(l1+1, l2+1)
		# Each position contains a 2-tuple of the respective gap score and True if the gap is a continuation, False if it is new
		self.leftMat = DirectionalMatrixWrapper(nrows=l1+1, ncols=l2+1) #numpy.zeros((l1+1, l2+1), dtype=('f16,b1')) 
		self.upMat = DirectionalMatrixWrapper(nrows=l1+1, ncols=l2+1) #numpy.zeros((l1+1, l2+1), dtype=('f16,b1'))
		self.scoreMat.set_data(0,0, 0)

		a_count = 0
		for i in range(1, l1 + 1): # set each row by the desired gap
			if self.composite != 2:
				# Update A(subtree) count; REMOVING - allowing gaps within a subtree for case of smaller tree matching part of larger tree
				#a_count += 1 if self.seq1.seq[i-1] == 'A' else -1 if self.seq1.seq[i-1] == 'T' else 0
				if a_count > 0: # Gap cannot start within a subtree
					self.scoreMat.set_data(i,0, None)
				else:
					self.scoreMat.set_data(i,0, self.costs['gap'] * i + self.costs['gapopen'])
			else:
				self.scoreMat.set_data(i,0, None)
			self.leftMat.score.set_data(i,0, None)
			self.upMat.score.set_data(i,0, None)
		a_count = 0
		for j in range(1, l2 + 1): # set each column by the desired gap
			if self.composite != 1:
				# Update A(subtree) count; REMOVING - allowing gaps within a subtree for case of smaller tree matching part of larger tree
				#a_count += 1 if self.seq2.seq[j-1] == 'A' else -1 if self.seq2.seq[j-1] == 'T' else 0
				if a_count > 0: # Gap cannot start within a subtree
					self.scoreMat.set_data(0,j, None)
				else:
					self.scoreMat.set_data(0,j, self.costs['gap'] * j + self.costs['gapopen'])
			else:
				self.scoreMat.set_data(0,j, None)
			self.leftMat.score.set_data(0,j, None)
			self.upMat.score.set_data(0,j, None)

		for i in range(1, l1+1): # per base-pair in sequence 1 ...
			for j in range(1, l2+1): # per base-pair in sequence 2, align them
			
				if (self.node_types[self.seq1.seq[i-1]] == 'C' and self.node_types[self.seq2.seq[j-1]] == 'A') or (self.node_types[self.seq1.seq[i-1]] == 'A' and self.node_types[self.seq2.seq[j-1]] == 'C'):
					score = None # no match if one is a C type and the other is an A type
				elif (self.node_types[self.seq1.seq[i-1]] == 'T') ^ (self.node_types[self.seq2.seq[j-1]] == 'T'):
					score = None # no match if one is a T type and the other is not
				elif self.scoreMat.get_data(i-1,j-1) is None:
					score = None # Diagonal is an unreachable position due to composite
				else:
					score = self.scoreMat.get_data(i-1,j-1) + get_score(self.seq1.seq[i-1], self.seq2.seq[j-1], self.submat)
				
				left, up = None, None
				if self.composite != 2: # If seq2 is composite, can't put gap characters in seq1
					# Cost for gapping left (over sequence 1)
					left, lefti, leftj = self.calculate_gap(i,
													j,
													self.seq1.seq,
													self.seq2.seq,
													self.scoreMat,
													self.directionMat,
													self.leftMat,
													self.TADict1,
													1)
				if self.composite != 1: # If seq1 is composite, can't put gap characters in seq2
					# Cost for gapping up (over sequence 2)
					up, upj, upi = self.calculate_gap(j,
												i,
												self.seq2.seq,
												self.seq1.seq,
												self.scoreMat.T,
												self.directionMat.T,
												self.upMat.T,
												self.TADict2,
												2)
				
				
				#stdout(str(i)+' '+str(j)+' match='+str(score)+' left='+str(left)+' up='+str(up))
				if score is not None and (left is None or score >= left) and (up is None or score >= up):
					# Node match is allowed and produces the best score
					self.scoreMat.set_data(i,j, score)
					self.directionMat.set_data(i,j, 0)
					self.backPos[i,j] = i-1,j-1
					#stdout('MATCH')
				elif left is not None and (up is None or left >= up):
					# Gapping left is allowed and produces the best score
					self.scoreMat.set_data(i,j, left)
					self.directionMat.set_data(i,j, 1)
					self.backPos[i,j] = lefti,leftj
					#stdout('LEFT')
				elif up is not None:
					# Gapping up is allowed and produces the best score
					self.scoreMat.set_data(i,j, up)
					self.directionMat.set_data(i,j, 2)
					self.backPos[i,j] = upi,upj
					#stdout('UP')
				else:
					# This location is unreachable due to presence of a composite sequence
					self.scoreMat.set_data(i,j,None)

		#### BACKTRACE ####
		i, j = l1, l2
		keepGapping = 0
		while i > 0 and j > 0: # walk-back to the index [0][0] of the m
			# Use next_i and next_j to keep track of backtrace position while holding on to current round's starting point
			next_i, next_j = i, j
			# if score is a gap in sequence 2 (direction is 1), only walk back on i
			if keepGapping == 1 or keepGapping == 0 and self.directionMat.get_data(i,j) == 1:
				# If the node being gapped is a T-node, appropriately gap the entire subtree
				if self.node_types[self.seq1.seq[i-1]] == 'T':
					# Get i_target and j_target depending on whether "forced" gapping or not
					if keepGapping == 1:
						# Determine whether T-gapping has includes the A node or matches A-C
						i_target = self.TADict1[i-1]
						if self.leftMat.a_c_match.get_data(i,j) == True:
							j_target = j-1
						else:
							j_target = j
					else:
						i_target = self.backPos[i,j][0]
						j_target = self.backPos[i,j][1]

					while next_i > i_target+1: # Walk back with gaps until the one position greater than the final position
						self.align1 += self.seq1.seq[next_i-1]
						self.align2 += '-'
						next_i -= 1

					# If j_target < next_j, then the gap is preceeded by an A-C match, so add that to the alignment
					if j_target < next_j:
						self.align1 += self.seq1.seq[next_i-1]
						self.align2 += self.seq2.seq[next_j-1]
						next_i -= 1
						next_j -= 1
					else: # otherwise the A is also gapped
						self.align1 += self.seq1.seq[next_i-1]
						self.align2 += '-'
						next_i -= 1

				# otherwise just gap the current node (it will be a C-node)
				else:
					self.align1 += self.seq1.seq[next_i-1]
					self.align2 += '-'
					next_i -= 1
					
				# Check whether the choice of gapping left requires the leftward position to also gap left
				if self.leftMat.extend_flag.get_data(i,j) == True:#[1]:
					keepGapping = 1
				else:
					keepGapping = 0
				
			# if score is a gap in sequence 1 (direction is 2), only walk back on j
			elif keepGapping == 2 or keepGapping == 0 and self.directionMat.get_data(i,j) == 2:
				if self.node_types[self.seq2.seq[j-1]] == 'T':
					# Get i_target and j_target depending on whether "forced" gapping or not
					if keepGapping == 2:
						# Determine whether T-gapping has includes the A node or matches A-C
						j_target = self.TADict2[j-1]
						if self.upMat.a_c_match.get_data(i,j) == True:
							i_target = i-1
						else:
							i_target = i
					else:
						i_target = self.backPos[i,j][0]
						j_target = self.backPos[i,j][1]
						
					while next_j > j_target+1:
						self.align1 += '-'
						self.align2 += self.seq2.seq[next_j-1]
						next_j -= 1
					if i_target < next_i:
						self.align1 += self.seq1.seq[next_i-1]
						self.align2 += self.seq2.seq[next_j-1]
						next_i -= 1
						next_j -= 1
					else:
						self.align1 += '-'
						self.align2 += self.seq2.seq[next_j-1]
						next_j -= 1
				else:
					self.align1 += '-'
					self.align2 += self.seq2.seq[next_j-1]
					next_j -= 1
					
				# Check whether the choice of gapping up requires the leftward position to also gap up
				if self.upMat.extend_flag.get_data(i,j) == True:
					keepGapping = 2
				else:
					keepGapping = 0

			# if the score is a match, walk-back one index in both i and j
			elif self.directionMat.get_data(i,j) == 0:
				keepGapping = 0
				self.align1 += self.seq1.seq[i-1]
				self.align2 += self.seq2.seq[j-1]
				next_i -= 1
				next_j -= 1
				
			i, j = next_i, next_j
		
		# walk-back to index 0 for both i and j; either could be reached first
		while i > 0:
			self.align1 += self.seq1.seq[i-1]
			self.align2 += '-'
			i -= 1
		while j > 0:
			self.align1 += '-'
			self.align2 += self.seq2.seq[j-1]
			j -= 1
		# Reverse the alignment strings as they are assembled backwards
		self.align1 = self.align1[::-1]
		self.align2 = self.align2[::-1]

# Global alignment of sequence to position weighted matrix, assuming pwm contains sequence
# This implementation only accepts A,C,T encoding and could be generalized
class PositionWeightedMatcher():
	def __init__(self, sequence, pwm, costs, node_types, allow_pwm_gaps=False, use_total=True):
		self.seq1 = sequence
		#self.seq2 = NeuriteSequence('PWM',pwm.node_type_sequence)
		self.seq2 = pwm.node_type_sequence
		self.costs = costs # dictionary of all costs (i.e. penalties)
		self.pwm = pwm.pwm # position weighted matrix / position specific score matrix: array of {char:weight} dictionary
		self.use_total = use_total
		self.create_node_types(node_types)
		self.TADict1 = None
		self.allow_pwm_gaps = allow_pwm_gaps
		if allow_pwm_gaps:
			self.TADict1 = create_ta_dictionary(self.seq1.seq, node_types, gap_cost=costs['gap'])
		self.TADict2 = create_ta_dictionary(self.seq2.seq, node_types, gap_cost=costs['gap'])
		if use_total:
			self.get_score = self.get_weight
		else:
			self.get_score = self.get_char_weight
		self.scoreMat = None # references score matrix
		self.directionMat = None # references diag(0),left(1),up(2) matrix
		self.upMat = None # references diag(0),left(1),up(2) matrix
		self.leftMat = None # references diag(0),left(1),up(2) matrix
		self.backPos = {} # the backtrace position from one position to its prior (contains integer pairs)
		self.align = '' # alignment string
		self.pwm_align = '' # pwm alignment string will be the composite with gaps if any
		self._aligner()
	
	def create_node_types(self,node_types):
		self.node_types = {}
		for node_type in node_types.keys():
			for residue in node_types[node_type]:
				self.node_types[residue] = node_type

	def get_char_weight(self,position, char1):
		return self.pwm[position][char1]

	# char1 is ignored and only used so that it takes the same parameters as get_char_weight
	def get_weight(self,position, char1='total'):
		return self.pwm[position]['total']
		
	# return top (highest) alignment score given sequence 1 and 2
	def get_top_score(self):
		return self.scoreMat.get_data(-1,-1)
	
	# alignment string resultant from sequence 1 and 2
	def get_alignment(self):
		return (self.align, self.pwm_align)
	
	# Create parser-friendly output given a NW alignment 
	def prettify(self):
		return [ self.get_top_score(), self.get_alignment(), self.seq2.name ]

	def determine_open_extend(self,i,j,m,directionM,dirScoreM,currentGapCost,gap_direction):
		if m.get_data(i,j) is None:
			return None,False
		gapScore = m.get_data(i,j) + currentGapCost
		# Add on the gap open cost if the prior position is not gapped in the same direction (gapping the same sequence)
		if dirScoreM.score.get_data(i,j) is None:#[0]): # Previous position can't gap
			gapScore += self.costs['gapopen']
			scoreExtendPair = gapScore,False
		elif directionM.get_data(i,j) is not gap_direction: # Previous position didn't choose gap in the same direction
			# Determine what the score would be if the previous position did gap in the same direction
			extendGapScore = dirScoreM.score.get_data(i,j) + currentGapCost
			openGapScore = gapScore + self.costs['gapopen']
			if openGapScore >= extendGapScore: # new gap is best, go with previous position's choice
				scoreExtendPair = openGapScore,False
			else: 
				# continuation is best, override previous position choice; if current position is on path, previous position will gap
				scoreExtendPair = extendGapScore,True
		else: # Previous position did choose gap, so this is automatically a continuation
			scoreExtendPair = gapScore,True
		return scoreExtendPair
		
	# Calculates the gap cost in a given direction from a given position, which depends on the node type
	def calculate_gap(self,i,j,seq1,seq2,m,directionM,dirScoreM,TADict,gap_direction):
		isExtend = False
		a_c_match = False
		# CType: gap one
		if self.node_types[seq1[i-1]] == 'C':
			# The prior position assuming a gap (index based on m)
			gapPosi = i - 1
			gapPosj = j
			# Determine the appropriate gap score depending on whether this opens or extends a gap
			gapScore,isExtend = self.determine_open_extend(gapPosi,
														gapPosj,
														m,
														directionM,
														dirScoreM,
														self.costs['gap'],
														gap_direction)

		# T-type: gap until paired A
		elif self.node_types[seq1[i-1]] == 'T': 
			# The prior position assuming a gap (index based on m)
			gapPosi = TADict[i-1]
			# Case where this is the last T; handle sentinal and get cost of front-gap
			if TADict[i-1] is -1:
				gapPosi = 0
			gapCostMajor = TADict[str(i-1)] # Cost of the gap from the T-node up to the A-node
			gapCostStart = self.costs['gap'] # Cost of the A-node that starts the gap

			# Calculate the total gap cost assuming the associated A is also gapped
			gapScoreGapFinish,isExtend = self.determine_open_extend(gapPosi,
																j,
																m,
																directionM,
																dirScoreM,
																gapCostMajor+gapCostStart,
																gap_direction)

			# If seq2 character is C-type, determine whether to match T-paired A and C, or to just gap the A
			# This will not happen if this is the last T (TADict[i-1] is -1), as the whole sequence must be gapped
			if self.node_types[seq2[j-1]] == 'C' and TADict[i-1] is not -1:
				# Calculate the total gap cost assuming the associated A-node matches a C-node
				if gap_direction == 1:
					ACScore = self.get_score(j-1,seq1[gapPosi])
				else:
					ACScore = self.get_score(gapPosi,seq2[j-1])

				if m.get_data(gapPosi,j-1) is None:
					gapScoreACFinish = None
				else:
					gapScoreACFinish = m.get_data(gapPosi,j-1) + gapCostMajor + ACScore + self.costs['gapopen']

				# Determine which gap produces a higher overall score, and use that for this position's gap score
				if gapScoreACFinish is not None and (gapScoreGapFinish is None or gapScoreACFinish >= gapScoreGapFinish):
					gapScore = gapScoreACFinish
					gapPosj = j-1
					isExtend = False
					a_c_match = True
				elif gapScoreGapFinish is not None:
					gapScore = gapScoreGapFinish
					gapPosj = j
					# isExtend depends on result of determine_open_extend
				else:
					gapScore = None
					gapPosj = j
					isExtend = False
			# If seq2 character is not a C, then use the total gap score assuming the A is also gapped
			else:
				gapScore = gapScoreGapFinish
				gapPosj = j
				# isExtend depends on result of determine_open_extend
		else: # A-type, no gapping allowed
			gapScore = None # original was set to None
			gapPosi = None
			gapPosj = None 
		dirScoreM.score.set_data(i,j, gapScore)
		dirScoreM.extend_flag.set_data(i,j, isExtend)
		dirScoreM.a_c_match.set_data(i,j, a_c_match)
		return [gapScore, gapPosi, gapPosj]

	# Execute alignment
	def _aligner(self):
		l1, l2 = len(self.seq1.seq), len(self.seq2.seq)	
		self.scoreMat = StateMatrix(l1+1, l2+1) # create matrix for storing counts
		self.directionMat = StateMatrix(l1+1, l2+1)
		# Each position contains a 2-tuple of the respective gap score and True if the gap is a continuation, False if it is new
		self.leftMat = DirectionalMatrixWrapper(nrows=l1+1, ncols=l2+1)
		self.upMat = DirectionalMatrixWrapper(nrows=l1+1, ncols=l2+1)
		self.scoreMat.set_data(0,0, 0)
		a_count = 0
		for i in range(1, l1 + 1): # set each row by the desired gap
			if self.allow_pwm_gaps:
				# Update A(subtree) count
				a_count += 1 if self.seq1.seq[i-1] == 'A' else -1 if self.seq1.seq[i-1] == 'T' else 0
				if a_count > 0: # Gap cannot start within a subtree
					self.scoreMat.set_data(i,0, None)
				else:
					self.scoreMat.set_data(i,0, self.costs['gap'] * i + self.costs['gapopen'])
			else:
				self.scoreMat.set_data(i,0, None)
			self.leftMat.score.set_data(i,0, None)
			self.upMat.score.set_data(i,0, None)
		a_count = 0
		for j in range(1, l2 + 1): # set each column by the desired gap
			# Update A(subtree) count
			a_count += 1 if self.seq2.seq[j-1] == 'A' else -1 if self.seq2.seq[j-1] == 'T' else 0
			if a_count > 0: # Gap cannot start within a subtree
				self.scoreMat.set_data(0,j, None)
			else:
				self.scoreMat.set_data(0,j, self.costs['gap'] * j + self.costs['gapopen'])

			#self.scoreMat.set_data(0,j, self.costs['gap'] * j + self.costs['gapopen'])
			self.leftMat.score.set_data(0,j, None)
			self.upMat.score.set_data(0,j, None)

		for i in range(1, l1+1): # per base-pair in sequence 1 ...
			for j in range(1, l2+1): # per base-pair in sequence 2, align them
				#print(str(i)+" "+str(j)+ " of "+str(l1)+" " +str(l2))
				if (self.node_types[self.seq1.seq[i-1]] == 'C' and self.node_types[self.seq2.seq[j-1]] == 'A') or (self.node_types[self.seq1.seq[i-1]] == 'A' and self.node_types[self.seq2.seq[j-1]] == 'C'):
					score = None # no match if one is a C type and the other is an A type
				elif (self.node_types[self.seq1.seq[i-1]] == 'T') ^ (self.node_types[self.seq2.seq[j-1]] == 'T'):
					score = None # no match if one is a T type and the other is not
				elif self.scoreMat.get_data(i-1,j-1) is None:
					score = None # Diagonal is an unreachable position due to composite
				else:
					score = self.scoreMat.get_data(i-1,j-1) + self.get_score(j-1,self.seq1.seq[i-1])

				left, up = None, None
				if self.allow_pwm_gaps: # Only calculate gap costs in this direction if gaps in pwm are allowed
					# Cost for gapping left (over sequence 1)
					left, lefti, leftj = self.calculate_gap(i,
										j,
										self.seq1.seq,
										self.seq2.seq,
										self.scoreMat,
										self.directionMat,
										self.leftMat,
										self.TADict1,
										1)

				# Cost for gapping up (over sequence 2)
				up, upj, upi = self.calculate_gap(j,
								  i,
								  self.seq2.seq,
								  self.seq1.seq,
								  self.scoreMat.T,
								  self.directionMat.T,
								  self.upMat.T,
								  self.TADict2,
								  2)
				
				
				#stdout(str(i)+' '+str(j)+' match='+str(score)+' left='+str(left)+' up='+str(up))
				if score is not None and (left is None or score >= left) and (up is None or score >= up):
					# Node match is allowed and produces the best score
					self.scoreMat.set_data(i,j, score)
					self.directionMat.set_data(i,j, 0)
					self.backPos[i,j] = i-1,j-1
					#stdout('MATCH')
				elif left is not None and (up is None or left >= up):
					# Gapping left is allowed and produces the best score
					self.scoreMat.set_data(i,j, left)
					self.directionMat.set_data(i,j, 1)
					self.backPos[i,j] = lefti,leftj
					#stdout('LEFT')
				elif up is not None:
					# Gapping up is allowed and produces the best score
					self.scoreMat.set_data(i,j, up)
					self.directionMat.set_data(i,j, 2)
					self.backPos[i,j] = upi,upj
					#stdout('UP')
				else:
					# This location is unreachable due to presence of a pre-consensus sequence
					self.scoreMat.set_data(i,j,None)


		i, j = l1, l2 # for trace-back process 
		keepGapping = 0
		while i > 0 and j > 0: # walk-back to the index [0][0] of the m
			# Use next_i and next_j to keep track of backtrace position while holding on to current round's starting point
			next_i, next_j = i, j

			#if j == 1106:
			#	print(self.seq1.seq[i-1]+'vs'+self.seq2.seq[j-1]+'; keepGapping: '+str(keepGapping)+'; acMatch: '+str(self.upMat.a_c_match.get_data(i,j)))

			# if score is a gap in sequence 2 (direction is 1), only walk back on i
			if keepGapping == 1 or keepGapping == 0 and self.directionMat.get_data(i,j) == 1:
				# If the node being gapped is a T-node, appropriately gap the entire subtree
				if self.node_types[self.seq1.seq[i-1]] == 'T':
					# Get i_target and j_target depending on whether "forced" gapping or not
					if keepGapping == 1:
						i_target = self.TADict1[i-1]
						if self.leftMat.a_c_match.get_data(i,j) == True:
							j_target = j-1
						else:
							j_target = j
					else:
						i_target = self.backPos[i,j][0]
						j_target = self.backPos[i,j][1]

					while next_i > i_target+1: # Walk back with gaps until the one position greater than the final position
						self.align += self.seq1.seq[next_i-1]
						self.pwm_align += '-'
						next_i -= 1
					if j_target < j: # If j_target < j, then the gap is preceeded by an A-C match, so add that to the alignment
						self.align += self.seq1.seq[next_i-1]
						self.pwm_align += self.seq2.seq[next_j-1]
						next_i -= 1
						next_j -= 1
					else: # otherwise the A is also gapped
						self.align += self.seq1.seq[next_i-1]
						self.pwm_align += '-'
						next_i -= 1
				# otherwise just gap the current node (it will be a C-node)
				else:
					self.align += self.seq1.seq[next_i-1]
					self.pwm_align += '-'
					next_i -= 1
					
				# Check whether the choice of gapping left requires the leftward position to also gap left
				if self.leftMat.extend_flag.get_data(i,j) == True:#[1]:
					keepGapping = 1
				else:
					keepGapping = 0
				
			# if score is a gap in sequence 1 (direction is 2), only walk back on j
			elif keepGapping == 2 or keepGapping == 0 and self.directionMat.get_data(i,j) == 2:
				if self.node_types[self.seq2.seq[j-1]] == 'T':
					# Get i_target and j_target depending on whether "forced" gapping or not
					if keepGapping == 2:
						if self.upMat.a_c_match.get_data(i,j) == True:
							i_target = i-1
						else:
							i_target = i
						j_target = self.TADict2[j-1]
					else:
						i_target = self.backPos[i,j][0]
						j_target = self.backPos[i,j][1]
						
					while next_j > j_target+1:
						self.align += '-'
						self.pwm_align += self.seq2.seq[next_j-1]
						next_j -= 1
					if i_target < next_i:
						self.align += self.seq1.seq[next_i-1]
						self.pwm_align += self.seq2.seq[next_j-1]
						next_i -= 1
						next_j -= 1
					else:
						self.align += '-'
						self.pwm_align += self.seq2.seq[next_j-1]
						next_j -= 1
				else:
					self.align += '-'
					self.pwm_align += self.seq2.seq[next_j-1]
					next_j -= 1
					
				# Check whether the choice of gapping up requires the leftward position to also gap up
				if self.upMat.extend_flag.get_data(i,j) == True:
					keepGapping = 2
				else:
					keepGapping = 0

			# if the score is a match, walk-back one index in both i and j
			elif self.directionMat.get_data(i,j) == 0:
				keepGapping = 0
				self.align += self.seq1.seq[i-1]
				self.pwm_align += self.seq2.seq[j-1]
				next_i -= 1
				next_j -= 1

			i, j = next_i, next_j
		
		# walk-back to index 0 for both i and j; either could be reached first
		while i > 0:
			self.align += self.seq1.seq[i-1]
			self.pwm_align += '-'
			i -= 1
		while j > 0:
			self.align += '-'
			self.pwm_align += self.seq2.seq[j-1]
			j -= 1
		# Reverse the alignment strings as they are assembled backwards
		self.align = self.align[::-1]
		self.pwm_align = self.pwm_align[::-1]

class LocalAlignmentWrapper():
	def __init__(self, s1, s2, align1, align2, score, start1=None, start2=None, end1=None, end2=None):
		self.s1 = s1
		self.s2 = s2
		self.align1 = align1
		self.align2 = align2
		self.score = score
		self.start1 = start1
		self.start2 = start2
		self.end1 = end1
		self.end2 = end2
		
# # Implementation of local alignment - Smith-Waterman
# class TreewiseSmithWaterman():
# 	def __init__(self, s1, s2, costs, submat, node_types):
# 		self.seq1 = s1 # sequence 1
# 		self.seq2 = s2 # sequence 2
# 		self.costs = costs # dictionary of all costs (i.e. penalties)
# 		self.submat = submat # substitution matrix
# 		self.create_node_types(node_types)
# 		self.create_residue_specific_gapcost()
# 		self.TADict1 = create_ta_dictionary(s1.seq, node_types, submat, costs['gap'])
# 		self.TADict2 = create_ta_dictionary(s2.seq, node_types, submat, costs['gap'])
# 		self.scoreMat = None # references score matrix
# 		self.directionMat = None # references diag(0),left(1),up(2) matrix
# 		self.leftMat = None # references diag(0),left(1),up(2) matrix
# 		self.upMat = None # references diag(0),left(1),up(2) matrix
# 		self.backPos = {} # the backtrace position from one position to its prior (contains integer pairs)
# 		self._alignments = []
# 		self.masked = None
# 		self.min_score = None
# 		self.curr_min_score = None
# 		self._aligner()
# 
# 	def create_node_types(self,node_types):
# 		self.node_types = {}
# 		for node_type in node_types.keys():
# 			for residue in node_types[node_type]:
# 				self.node_types[residue] = node_type
# 	
# 	# Fills in all gap-residue pairs either with flipped order entry if it exists, else with the default gap cost
# 	# Also fills in flipped residue-residue scores
# 	def create_residue_specific_gapcost(self):
# 		newToSubmat = {}
# 		for pair in self.submat.keys():
# 			if pair[0] == '-' and pair[1] == '-':
# 				# do nothing, this is useless and shouldn't happen
# 				pass
# 			elif pair[0] == '-' or pair[1] == '-':
# 				# Note that the given residue has an associated gap cost
# 				if pair[0] == '-' and (pair[1],'-') not in self.submat.keys():
# 					newToSubmat[pair[1],'-'] = self.submat[pair]
# 				elif pair[1] == '-' and (pair[1],'-') not in self.submat.keys():
# 					newToSubmat[pair[1],'-'] = self.submat[pair]
# 			else:
# 				# Fill the residueDict so none are missed
# 				if (pair[0],'-') not in self.submat.keys() and ('-',pair[0]) not in self.submat.keys():
# 					newToSubmat[pair[0],'-'] = self.costs['gap']
# 					newToSubmat['-',pair[0]] = self.costs['gap']
# 				if (pair[1],'-') not in self.submat.keys() and ('-',pair[1]) not in self.submat.keys():
# 					newToSubmat[pair[1],'-'] = self.costs['gap']
# 					newToSubmat['-',pair[1]] = self.costs['gap']
# 				if (pair[1],pair[0]) not in self.submat.keys():
# 					newToSubmat[pair[1],pair[0]] = self.submat[pair]
# 		for pair in newToSubmat.keys():
# 			self.submat[pair] = newToSubmat[pair]
# 	
# 	def determine_open_extend(self,i,j,m,directionM,dirScoreM,currentGapCost,gap_direction):
# 		if m.get_data(i,j) is None:
# 			return None,False
# 		gapScore = m.get_data(i,j) + currentGapCost
# 		# Add on the gap open cost if the prior position is not gapped in the same direction (gapping the same sequence)
# 		if dirScoreM.score.get_data(i,j) is None:#[0]): # Previous position can't gap
# 			gapScore += self.costs['gapopen']
# 			scoreExtendPair = gapScore,False
# 		elif directionM.get_data(i,j) is not gap_direction: # Previous position didn't choose gap
# 			extendGapScore = dirScoreM.score.get_data(i,j) + currentGapCost
# 			openGapScore = gapScore + self.costs['gapopen']
# 			if openGapScore >= extendGapScore: # new gap is best, go with previous position's choice
# 				scoreExtendPair = openGapScore,False
# 			else: 
# 				# continuation is best, override previous position choice; if current position is on path, previous position will gap
# 				scoreExtendPair = extendGapScore,True
# 		else: # Previous position did choose gap, so this is automatically a continuation
# 			scoreExtendPair = gapScore,True
# 		return scoreExtendPair
# 		
# 	# Calculates the gap cost in a given direction from a given position, which depends on the node type
# 	def calculate_gap(self,i,j,seq1,seq2,m,directionM,dirScoreM,TADict,gap_direction):
# 		isExtend = False
#		a_c_match = False
# 		# CType: gap one
# 		if self.node_types[seq1[i-1]] == 'C':
# 			# The prior position assuming a gap (index based on m)
# 			gapPosi = i - 1
# 			gapPosj = j
# 			# Determine the appropriate gap score depending on whether this opens or extends a gap
# 			gapScore,isExtend = self.determine_open_extend(gapPosi,
# 														gapPosj,
# 														m,
# 														directionM,
# 														dirScoreM,
# 														get_gapcost(seq1[i-1],self.submat),
# 														gap_direction)
# 
# 		# TType: gap until paired A
# 		elif self.node_types[seq1[i-1]] == 'T': 
# 			# The prior position assuming a gap (index based on m)
# 			gapPosi = TADict[i-1]
# 			# Case where this is the last T; handle sentinal and get cost of front-gap
# 			if TADict[i-1] is -1:
# 				gapPosi = 0
# 			gapCostMajor = TADict[str(i-1)] # Cost of the gap from the T-node up to the A-node
# 			gapCostStart = get_gapcost(seq1[gapPosi],self.submat) # Cost of the A-node that starts the gap
# 
# 			# Calculate the total gap cost assuming the associated A is also gapped
# 			gapScoreGapFinish,isExtend = self.determine_open_extend(gapPosi,
# 																j,
# 																m,
# 																directionM,
# 																dirScoreM,
# 																gapCostMajor+gapCostStart,
# 																gap_direction)
# 
# 			# If seq2 character is C-type, determine whether to match T-paired A and C, or to just gap the A
# 			# This will not happen if this is the last T (TADict[i-1] is -1), as the whole sequence must be gapped
# 			if self.node_types[seq2[j-1]] == 'C' and TADict[i-1] is not -1:
# 				# Calculate the total gap cost assuming the associated A-node matches a C-node
# 				ACScore = get_score(seq1[gapPosi],seq2[j-1],self.submat)
# 				if m.get_data(gapPosi,j-1) is None:
# 					gapScoreACFinish = None
# 				else:
# 					gapScoreACFinish = m.get_data(gapPosi,j-1) + gapCostMajor + ACScore + self.costs['gapopen']
# 
# 				# Determine which gap produces a higher overall score, and use that for this position's gap score
# 				if gapScoreACFinish is not None and (gapScoreGapFinish is None or gapScoreACFinish >= gapScoreGapFinish):
# 					gapScore = gapScoreACFinish
# 					gapPosj = j-1
# 					isExtend = False
#					a_c_match = True
# 				elif gapScoreGapFinish is not None:
# 					gapScore = gapScoreGapFinish
#					gapPosj = j
# 				else:
# 					gapScore = None
# 					gapPosj = j
# 					isExtend = False
# 			# If seq2 character is not a C, then use the total gap score assuming the A is also gapped
# 			else:
# 				gapScore = gapScoreGapFinish
# 				gapPosj = j
# 		else: # AType, no gapping allowed
# 			gapScore = None # original was set to None
# 			gapPosi = None
# 			gapPosj = None 
# 		dirScoreM.score.set_data(i,j, gapScore)
# 		dirScoreM.extend_flag.set_data(i,j, isExtend)
#		dirScoreM.a_c_match.set_data(i,j, a_c_match)
# 		return [gapScore, gapPosi, gapPosj]
# 
# 	# Execute alignment
# 	def _aligner(self):
# 		l1, l2 = len(self.seq1.seq), len(self.seq2.seq)	
# 		self.scoreMat = Matrix(l1+1, l2+1) # create matrix for storing counts
# 		self.directionMat = Matrix(l1+1, l2+1)
# 		# Each position contains a 2-tuple of the respective gap score and True if the gap is a continuation, False if it is new
# 		self.leftMat = DirectionalMatrixWrapper(nrows=l1+1, ncols=l2+1) #numpy.zeros((l1+1, l2+1), dtype=('f16,b1')) 
# 		self.upMat = DirectionalMatrixWrapper(nrows=l1+1, ncols=l2+1) #numpy.zeros((l1+1, l2+1), dtype=('f16,b1'))
# 		self.scoreMat.set_data(0,0, 0)
# 		for i in range(1, l1 + 1): # set each row by the desired gap
# 			self.scoreMat.set_data(i,0, 0)
# 			self.leftMat.score.set_data(i,0, None) #[0] = None
# 			self.upMat.score.set_data(i,0, None)
# 		for j in range(1, l2 + 1): # set each column by the desired gap
# 			self.scoreMat.set_data(0,j, 0)
# 			self.leftMat.score.set_data(0,j, None) #[0] = None
# 			self.upMat.score.set_data(0,j, None)
# 			
# 		self.masked = DirectionalMatrixWrapper(nrows=l1+1, ncols=l2+1)
# 
# 		for i in range(1, l1+1): # per base-pair in sequence 1 ...
# 			for j in range(1, l2+1): # per base-pair in sequence 2, align them
# 				self.mask.set_data(i,j,1)
# 				if (self.node_types[self.seq1.seq[i-1]] == 'C' and self.node_types[self.seq2.seq[j-1]] == 'A') or (self.node_types[self.seq1.seq[i-1]] == 'A' and self.node_types[self.seq2.seq[j-1]] == 'C'):
# 					score = None # no match if one is a C type and the other is an A type
# 				elif (self.node_types[self.seq1.seq[i-1]] == 'T') ^ (self.node_types[self.seq2.seq[j-1]] == 'T'):
# 					score = None # no match if one is a T type and the other is not
# 				elif self.scoreMat.get_data(i-1,j-1) is None:
# 					score = None # Diagonal is an unreachable position due to pre-consensus
# 				else:
# 					score = self.scoreMat.get_data(i-1,j-1) + get_score(self.seq1.seq[i-1], self.seq2.seq[j-1], self.submat)
# 				
# 				# Cost for gapping left (over sequence 1)
# 				left, lefti, leftj = self.calculate_gap(i,
# 											j,
# 											self.seq1.seq,
# 											self.seq2.seq,
# 											self.scoreMat,
# 											self.directionMat,
# 											self.leftMat,
# 											self.TADict1,
# 											1)
# 
# 				# Cost for gapping up (over sequence 2)
# 				up, upj, upi = self.calculate_gap(j,
# 											i,
# 											self.seq2.seq,
# 											self.seq1.seq,
# 											self.scoreMat.T,
# 											self.directionMat.T,
# 											self.upMat.T,
# 											self.TADict2,
# 											2)
# 				
# 			#	if (i == 539 and j == 31):
# 			#		print(str(i)+'; '+str(j)+ '; 1: ' +self.seq1.seq[i-1]+ '; 2: '+self.seq2.seq[j-1]+'; Sc: ' +str(score)+ '; Left: '+str(left)+'; Up: '+str(up))
# 			#	if (i == 538 and j == 30):
# 			#		print(str(i)+'; '+str(j)+ '; 1: ' +self.seq1.seq[i-1]+ '; 2: '+self.seq2.seq[j-1]+'; Sc: ' +str(score)+ '; Left: '+str(left)+'; Up: '+str(up))
# 				
# 				#stdout(str(i)+' '+str(j)+' match='+str(score)+' left='+str(left)+' up='+str(up))
# 				if score is not None and (left is None or score >= left) and (up is None or score >= up):
# 					# Node match is allowed and produces the best score
# 					self.scoreMat.set_data(i,j, score)
# 					self.directionMat.set_data(i,j, 0)
# 					self.backPos[i,j] = i-1,j-1
# 					#if (i == 539 and j == 31):
# 					#	print("Picked diagonal; prev pos is: "+str(self.backPos[i,j]))
# 					#stdout('MATCH')
# 				elif left is not None and (up is None or left >= up):
# 					# Gapping left is allowed and produces the best score
# 					self.scoreMat.set_data(i,j, left)
# 					self.directionMat.set_data(i,j, 1)
# 					self.backPos[i,j] = lefti,leftj
# 					#stdout('LEFT')
# 				elif up is not None:
# 					# Gapping up is allowed and produces the best score
# 					self.scoreMat.set_data(i,j, up)
# 					self.directionMat.set_data(i,j, 2)
# 					self.backPos[i,j] = upi,upj
# 					#stdout('UP')
# 				else:
# 					# This location is unreachable due to presence of a pre-consensus sequence
# 					self.scoreMat.set_data(i,j,None)
# 
# 				# SmithWaterman minimum cell score is 0
# 				if self.scoreMat.get_data(i,j) is not None and self.scoreMat.get_data(i,j) < 0:
# 					self.scoreMat.set_data(i,j,0)
# 
# 	def _find_next_alignment(self,min_score):
# 		bestScore = min_score
# 		bestIJ = 0,0
# 		l1, l2 = len(self.seq1.seq), len(self.seq2.seq)	
# 		for i in range(1, l1+1):
# 			for j in range(1, l2+1):
# 				if not self.masked.get_data(i,j) and self.scoreMat.get_data(i,j) >= bestScore:
# 					bestScore = self.scoreMat.get_data(i,j)
# 					bestIJ = i,j
# 		return bestIJ
# 	
# 	def _backtrace(self,i,j):
# 		keepGapping = 0
# 		align1 = ""
# 		align2 = ""
# 		initialI = i
# 		initialJ = j
# 		# walk-back until score of 0 is reached
# 		while self.scoreMat.get_data(i,j) > 0: 
# 			i_target = i
# 			j_target = j
# 			# Mask this position from future alignment searching
# 			self.masked.set_data(i,j, 1)
#
#                       next_i, next_j = i, j 
#
# 			# if score is a gap in sequence 2 (direction is 1), only walk back on i
# 			if keepGapping == 1 or keepGapping == 0 and self.directionMat.get_data(i,j) == 1:
# 				
# 				# If the node being gapped is a T-node, appropriately gap the entire subtree
# 				if self.node_types[self.seq1.seq[i-1]] == 'T':
# 					# Get i_target and j_target depending on whether "forced" gapping or not
# 					if keepGapping == 1:
# 						i_target = self.TADict1[i-1]
#						if self.leftMat.a_c_match.get_data(i,j) == True:
#							j_target = j-1
#						else:
#							j_target = j
# 					else:
# 						i_target = self.backPos[i,j][0]
# 						j_target = self.backPos[i,j][1]
# 
# 					#print(str(self.backPos[i,j]) + ' ' + str(i_target)+' '+str(j_target))
# 					while next_i > i_target+1: # Walk back with gaps until the one position greater than the final position
# 						align1 += self.seq1.seq[next_i-1]
# 						align2 += '-'
# 						next_i -= 1
# 					if j_target < j: # If j_target < j, then the gap is preceeded by an A-C match, so add that to the alignment
# 						align1 += self.seq1.seq[next_i-1]
# 						align2 += self.seq2.seq[next_j-1]
# 						next_i -= 1
# 						next_j -= 1
# 					else: # otherwise the A is also gapped
# 						align1 += self.seq1.seq[next_i-1]
# 						align2 += '-'
# 						next_i -= 1
# 				# otherwise just gap the current node (it will be a C-node)
# 				else:
# 					align1 += self.seq1.seq[next_i-1]
# 					align2 += '-'
# 					next_i -= 1
# 					
# 				# Check whether the choice of gapping left requires the leftward position to also gap left
# 				if self.leftMat.extend_flag.get_data(i,j) == True:#[1]:
# 					keepGapping = 1
# 				else:
# 					keepGapping = 0
# 
# 			# if score is a gap in sequence 1 (direction is 2), only walk back on j
# 			elif keepGapping == 2 or keepGapping == 0 and self.directionMat.get_data(i,j) == 2:
# 				if self.node_types[self.seq2.seq[j-1]] == 'T':
# 					# Get i_target and j_target depending on whether "forced" gapping or not
# 					if keepGapping == 2:
#						if self.upMat.a_c_match.get_data(i,j) == True:
#							i_target = i-1
#						else:
#							i_target = i
# 						j_target = self.TADict2[j-1]
# 					else:
# 						i_target = self.backPos[i,j][0]
# 						j_target = self.backPos[i,j][1]
# 
# 					while next_j > j_target+1:
# 						align1 += '-'
# 						align2 += self.seq2.seq[next_j-1]
# 						next_j -= 1
# 					if i_target < next_i:
# 						align1 += self.seq1.seq[next_i-1]
# 						align2 += self.seq2.seq[next_j-1]
# 						next_i -= 1
# 						next_j -= 1
# 					else:
# 						align1 += '-'
# 						align2 += self.seq2.seq[next_j-1]
# 						next_j -= 1
# 				else:
# 					align1 += '-'
# 					align2 += self.seq2.seq[next_j-1]
# 					next_j -= 1
# 					
# 				# Check whether the choice of gapping up requires the leftward position to also gap up
# 				if self.upMat.extend_flag.get_data(i,j) == True:
# 					keepGapping = 2
# 				else:
# 					keepGapping = 0
# 
# 			# if the score is a match, walk-back one index in both i and j
# 			elif self.directionMat.get_data(i,j) == 0:
# 				keepGapping = 0
# 				align1 += self.seq1.seq[i-1]
# 				align2 += self.seq2.seq[j-1]
# 				next_i -= 1
# 				next_j -= 1
# 		
#                       i, j = next_i, next_j
#
# 		# Reverse the alignment strings as they are assembled backwards
# 		align1 = align1[::-1]
# 		align2 = align2[::-1]
# 		
# 		alignment = LocalAlignmentWrapper(self.seq1,self.seq2,align1,align2,self.scoreMat.get_data(initialI,initialJ),i_target,j_target,initialI,initialJ)
# 		return alignment
# 
# 	def get_alignments(self,min_score):
# 		if self.min_score is None or min_score < self.min_score:
# 			found_new = true
# 			# Run alignment search and backtraces
# 			while found_new:
# 				# Find current best score and end alignment position
# 				i,j = _find_next_alignment(min_score)
# 				
# 				if i > 0: # An alignment better than min_score was found
# 					# Traceback to get alignment
# 					alignment = _backtrace(i,j)
# 					_alignments.extend(alignment)
# 					found_new = True
# 				else:
# 					found_new = False
# 
# 		# Go through current alignments and pull stdout those with sufficient score to return
# 		select_alignments = []
# 		for alignment in self._alignments:
# 			if alignment.score >= min_score:
# 				select_alignments.extend(alignment)
# 				
# 		self.min_score = min_score
# 		return select_alignments

# Maps the aligned two bases against a user-selected substitution matrix
def get_score(cA, cB, submatrix):
	return submatrix[(cA, cB)]
#	if (cA, cB) in submatrix:
#		return submatrix[(cA, cB)]
#	else:
#		return submatrix[(cB, cA)] # returns score

# # Determines whether the given character pair is in the substitution matrix
# def has_score(cA, cB, submatrix):
# 	return (cA, cB) in submatrix or (cB, cA) in submatrix

# Returns the gap cost of the given character
def get_gapcost(char,submatrix):
	return get_score(char,'-',submatrix)

