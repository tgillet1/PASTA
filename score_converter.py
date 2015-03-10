import concurrent.futures, numpy
from datetime import datetime

# Maps each query sequence against a set of targets (itself)
def update_scores(input_line, targets, sequence_lengths, costs):
	results = [] # K => target, V => aligned queries 
	
	query,scores_str = input_line.split('\t',1)
	query_length = sequence_lengths[query]
	scores = scores_str.split()
	normalized_scores = numpy.zeros(len(scores))

	#print('NumScores: '+str(len(scores)))

	for index in range(0,len(scores)):
		target_length = sequence_lengths[targets[index].name]

		if scores[index] == 'None':
			normalized_scores[index] = None
		else:
			normalized_score = int(float(scores[index]))
			if (query_length != target_length):
				normalized_score -= (abs(query_length-target_length)*costs['gap'] + costs['gapopen'])
			normalized_scores[index] = normalized_score/min(query_length,target_length)
	return query, normalized_scores

def generate_querylength_dictionary(queries):
	query_dict = {}
	for query in queries:
		query_dict[query.name] = len(query.seq)
	return query_dict
	
# Executes the pairwise application
class ScoreConversionDriver():
	def __init__(self, input_state, queries, targets=None):
		self.queries = queries # factory operates given targets and queries
		self.sequence_lengths = generate_querylength_dictionary(queries)
		if targets is None:
			self.targets = queries
		else:
			self.targets = targets
			if not targets == queries:
				self.sequence_lengths.update(generate_querylength_dictionary(targets))

		self.input_handle = open(input_state.get_args()['rawscore_file'],'r')
		self.costs = input_state.get_penalties() # set costs to factory

		self.score_handle = open(input_state.get_args()['output'], 'w') # output file
		
		self.num_workers = input_state.get_args()['num_processes']
		self.num_complete = 0

	# Initialize the factory given query sequences and input arguments
	def start(self):
		executor = concurrent.futures.ProcessPoolExecutor(self.num_workers)
		
		#queries = self.input_handle.readLine().strip().split('\t')
		try:
			first_line = True
			for line in self.input_handle:
				if first_line:
					#self.queries = line.strip().split()
					self.score_handle.write(line)
					self.score_handle.flush()
					first_line = False
				else:
					f = executor.submit(update_scores, line, self.queries, self.sequence_lengths, self.costs)
					f.add_done_callback(self._callback)
			executor.shutdown()
			self.close_output_buffers()
			print('** Analysis Complete **')
		except KeyboardInterrupt:
			executor.shutdown()

	# Close all I/O buffers such as file handles
	def close_output_buffers(self):
		self.score_handle.close()
		self.input_handle.close()

	# Get the headers, i.e. top-most row for the score matrix
	def _create_header(self, results):
		h = '\t' +'\t'.join([h[-1] for h in results])
		self.score_handle.write(h + '\n')
		self.score_handle.flush()

	def _callback(self, return_val):
		cbexcept = return_val.exception()
		if cbexcept is not None:
			print("Err1: "+str(cbexcept))
                        #print("Targets: "+str(list([x.seq for x in self.targets])))
                        #print("Queries: "+str(list([x.seq for x in self.queries])))

		res = return_val.result() # get result once thread is complete
		target, results = res
		
		# save scores to the alignment matrix
		scores = '\t'.join(["%.4f" % round(s,4) for s in results])
		self.score_handle.write(target + '\t' + scores + '\n')
		self.score_handle.flush()

		self.num_complete += 1
		#out(' --> ' + target + ' [OK] '+datetime.time(datetime.now())) # print-out progress
		#print(' --> ' + target + ' [OK] '+str(self.num_complete)+' of '+str(len(self.queries))+' at '+str(datetime.time(datetime.now()))) # print-out progress
