import re
from matrix import ContingencyMatrix
from buffer import status_message
import concurrent.futures
from parameter import generate_sequence_set
from msa import MultipleSequenceDriver, ConsensusFilterFactory
import sequence
import warnings

class DomainSetBuilder():
    ''' 
    Given arguments relative to domain analysis, this class enables the ability
    to analyze a consensus sequence and extract out domains which satisfy
    user-provided arguments. If allowable_treeseq_types is not set (None), this
    means no type test needs to be run (all types are allowed).
    '''
    def __init__(self, consensus, win, max_gap, is_strip, is_enum=False, allowable_treeseq_types=None, node_types=sequence.default_nodetypes, min_win=1):
        if isinstance(consensus,str):
            self.consensus_seq = consensus
        else: # is NeuriteSequence
            self.consensus_seq = consensus.seq
        if is_strip:
            self.consensus_seq = self.consensus_seq.replace('-','')

        self.is_enumerate = is_enum # enumerate window size and max #/gaps
        self.win = win # sliding-window size
        self.min_win = min_win
        self.max_gap = max_gap # maximum #/gaps in each domain
        #self.is_strip = is_strip
        self.allowable_treeseq_types = allowable_treeseq_types
        self.node_types = node_types
        
    def build(self):
        wins = [self.win]
        if self.is_enumerate: # enumeration yields dynamic window and gaps
            wins = list(range(self.min_win, self.win + 1)) # window sizes from 1 .. win
            
        # iterate over window and gaps (if enumeration) and identify domains
        domains = {} # key => domain, value => count (abundance)
        for w in wins:
            if w > len(self.consensus_seq):
                if not self.is_enumerate:
                    raise IOError('-win must be less than consensus length ('+str(len(self.consensus))+')')
                elif w == max(wins):
                    warnings.warn('-win ('+str(self.win)+') is larger than consensus length ('+str(len(self.consensus_seq))+'), so only domains of consensus length or smaller will be found.')
            else: # for each window, pull-out the respective domain
                for idx in range(len(self.consensus_seq)):
                    sub_str = self.consensus_seq[idx: idx + w] # reference
                    num_gap = sub_str.count('-') # count number of gaps
                    # if required, candidate domain must be a complete tree
                    if self.allowable_treeseq_types is None or sequence.tree_sequence_type(sub_str,self.node_types) in self.allowable_treeseq_types:
                        # only-gapped sequences are ignored (if stripped, there will be no gaps)
                        if not re.match('^-+$', sub_str) and num_gap <= self.max_gap:
                            if sub_str not in domains:
                                domains[sub_str] = 0
                            domains[sub_str] += 1 # increment its abundance

                    #if self.is_strip: # if gaps are found, remove them
                    #    sub_str = sub_str.replace('-', '')
                    #num_gap = sub_str.count('-') # count number of gaps
                    # domain must equal sliding window and, if required, a complete tree
                    #if len(sub_str) == w and (not self.complete_tree_f or is_complete_tree(sub_str)):
                    #    # only-gapped sequences are ignored
                    #    if not re.match('^-+$', sub_str) and num_gap <= self.max_gap:
                    #        if sub_str not in domains:
                    #            domains[sub_str] = 0
                    #        domains[sub_str] += 1 # increment its abundance

        return domains # return dictionary of domains and their abundances

    def is_complete_tree(self, seq):
        aCount = 0
        for c in seq:
            if c == 'A':
                aCount+=1
            elif c == 'T':
                aCount-=1
        if aCount == -1:
            return True

    def is_partial_tree(self, seq):
        aCount = 0
        for c in seq:
            if c == 'A':
                aCount+=1
            elif c == 'T':
                aCount-=1
        if aCount == -1:
            return True


class DomainAbundanceBuilder():
    ''' 
    This class enables the ability to build contingency matrices so as to
    enable the ability to contrast abundance of an item across different 
    groups. This then leads to the fact that a domain, D, could be abundant in
    1 of 3 possible scenarios: 
    1) D is present in query, Q, and not baseline, B.
    2) D is present in baseline, B, and not query, Q.
    3) D is present in both baseline, B, and query, Q.
    Thus, for all scenarios, the contingency matrix must model group-specific
    abundances so that truly over-represented domains can be identified.
    '''
    def __init__(self, query, baseline):
        self.query = query # query domain abundances
        self.baseline = baseline # baseline domain abundances
    
    def _query_exclusive_domains(self):
        names = set()
        for q in self.query:
            if q not in self.baseline:
                names.add(q)
        return names # references domains only found in the query
    
    def _baseline_exclusive_domains(self):
        names = set()
        for b in self.baseline:
            if b not in self.query:
                names.add(b)
        return names # references domains only found in the baseline
    
    def _intersection_domains(self):
        # get intersection of domains present in both query and baseline sets
        return set(self.query.keys()).intersection(self.baseline.keys())
    
    def build(self):
        # Suppose we set the following: 
        # Query (G) and Baseline (! G)
        # Domain (i) and all-other domains (! i) 
        size_query = sum(self.query.values()) # same as n(G)
        size_baseline = sum(self.baseline.values()) # same as n(! G)
        matrices = []
        
        # building matrices in query and baseline
        for d in self._intersection_domains():
            i_and_G = self.query[d] # domain count in query (idx: 0, 0)
            i_and_not_G = self.baseline[d] # domain count in query (idx: 0, 1)
            not_i_and_G = size_query - i_and_G # not-domain in query (idx: 1, 0)
            not_i_and_not_G = size_baseline - i_and_not_G # not-domain in baseline (idx: 1, 1)
            cm = ContingencyMatrix(node = d, i_g = i_and_G, 
                                   i_not_g = i_and_not_G, not_i_g=not_i_and_G, 
                                   not_i_not_g = not_i_and_not_G)
            matrices.append(cm)
            
        # building matrices in baseline only
        for d in self._baseline_exclusive_domains():
            i_and_G = 0.01 # domain count in query (idx: 0, 0)
            i_and_not_G = self.baseline[d] # domain count in query (idx: 0, 1)
            not_i_and_G = size_query - i_and_G # not-domain in query (idx: 1, 0)
            not_i_and_not_G = size_baseline - i_and_not_G # not-domain in baseline (idx: 1, 1)
            cm = ContingencyMatrix(node = d, i_g = i_and_G, 
                                   i_not_g = i_and_not_G, not_i_g=not_i_and_G, 
                                   not_i_not_g = not_i_and_not_G)
            matrices.append(cm)
        
        # building matrices in query only
        for d in self._query_exclusive_domains():
            i_and_G = self.query[d] # domain count in query (idx: 0, 0)
            i_and_not_G = 0.01 # domain count in query (idx: 0, 1)
            not_i_and_G = size_query - i_and_G # not-domain in query (idx: 1, 0)
            not_i_and_not_G = size_baseline - i_and_not_G # not-domain in baseline (idx: 1, 1)
            cm = ContingencyMatrix(node = d, i_g = i_and_G, 
                                   i_not_g = i_and_not_G, not_i_g=not_i_and_G, 
                                   not_i_not_g = not_i_and_not_G)
            matrices.append(cm)
        
        return matrices
    
class DomainPrettyPrinter():
    ''' 
    Prints statistically-significant domains after computing a hypergeometric
    p-value representing the abundance for it.
    @param domains: List of identified domains
    @param pval: User-provided p-value cutoff
    '''
    def __init__(self, domains, pval, out):
        self.domains = domains
        self.pval = pval
        self.fname = out
    
    def display(self):
        ''' 
        Prints statistically-significant domains to the screen.
        '''
        handle = open(self.fname, 'w')
        handle.write('Domain\tp-value\n') # write header
        handle.flush()
        for i in self.domains:
            dom_pval = i.get_hypergeometric_pval() # compute domain p-value
            if dom_pval <= self.pval: # domain must be less than p-value cutoff
                handle.write(i.name + '\t' + str(round(dom_pval, 6)) + '\n')
                handle.flush()
        handle.close()
        

def _extract_domains(targets, is_baseline, input_state, threshold, thresh_type, max_domain_size=15, min_domain_size=1):
    # Run MSA
    msa_driver = MultipleSequenceDriver(targets, input_state)
    msa_driver.build_composite()
    msa_driver.align()

    # Derive consensus
    consensus_obj = msa_driver.build_consensus(threshold,thresh_type)
    
    #consensus_fact = ConsensusFilterFactory(msa_driver.alns,msa_driver.composite, threshold, thresh_type)
    #consensus_fact.build_consensus()
    #if (is_baseline):
    #    print('Baseline Consensus: '+str(consensus_fact.consensus).replace('-',''))
    #else:
    #    print('Target Consensus: '+str(consensus_fact.consensus).replace('-',''))

    # Extract domains (consensus_obj.consensus is already stripped of - chars)
    domainBuilder = DomainSetBuilder(consensus_obj.consensus,max_domain_size,0,True,is_enum=True,allowable_treeseq_types=['complete_tree','incomplete_tree'],min_win=min_domain_size)

    domains = domainBuilder.build()
    return is_baseline, domains, consensus_obj.consensus

def merge_counts(counts1, counts2):
    merged_counts = {}
    for key in counts1.keys():
        merged_counts[key] = counts1[key]
    for key in counts2.keys():
        if key in merged_counts.keys():
            merged_counts[key] += counts2[key]
        else:
            merged_counts[key] = counts2[key]
    return merged_counts

class DomainExtractionDriver():
    '''
    Runs domain extraction multiple times using either different sequence order and/or subset on both
    target and baseline sets, then compiles the list of unique domains and how frequently each occurs
    in either set.
    '''
    def __init__(self,targets,baselines,node_types,subsmat,num_runs,input_state,debug=0):
        self.targets = targets
        self.baselines = baselines
        self.node_types = node_types
        self.subsmat = subsmat
        self.num_runs = num_runs
        self.input_state = input_state
        self.args = input_state.get_args()
        self.domain_set = set() # domains from targets and baselines
        self.target_domains = {}
        self.baseline_domains = {}
        self.target_consensuses = []
        self.baseline_consensuses = []
        self.debug=debug
    
    def start(self):
        status_message('Multiple alignment running','please wait')
        args = self.args
        executor = concurrent.futures.ProcessPoolExecutor(args['n'])
        try:
            for i in range(self.num_runs):
                if self.debug >= 1:
                    print("DomainExtractionDriver: Run "+str(i+1)+" of "+str(self.num_runs))
                targets_sub = generate_sequence_set(self.targets,args['subsample'],args['random_subset'],args['random_order'],args['subsample_start'])
                f = executor.submit(_extract_domains, targets_sub, is_baseline=False,
                                    input_state=self.input_state, threshold=args['thresh'],
                                    thresh_type=args['type'], max_domain_size=args['win'], min_domain_size=args['minwin'])
                f.add_done_callback(self._callback)

                disallowed = []
                if args['disjoint_subset']:
                    disallowed = targets_sub
                baselines_sub = generate_sequence_set(self.baselines,args['subsample'],args['random_subset'],args['random_order'],args['subsample_start'],disallowed=disallowed)
                f = executor.submit(_extract_domains, baselines_sub, is_baseline=True,
                                    input_state=self.input_state, threshold=args['thresh'],
                                    thresh_type=args['type'], max_domain_size=args['win'], min_domain_size=args['minwin'])
                f.add_done_callback(self._callback)

            executor.shutdown()
            #self.close_output_buffers()
            status_message('Analysis complete', 'OK')
        except KeyboardInterrupt:
            executor.shutdown()

    def _callback(self, return_val):
        cbexcept = return_val.exception()
        if cbexcept is not None:
            print("Err1: "+str(cbexcept))

        is_baseline, domains, consensus = return_val.result()

        if self.debug >= 1:
            print("Completed a DED run")
        
        # Merge domain counts for current iteration with those from other iterations
        if is_baseline:
            self.baseline_domains = merge_counts(self.baseline_domains,domains)
            self.baseline_consensuses.append(consensus)
        else:
            self.target_domains = merge_counts(self.target_domains,domains)
            self.target_consensuses.append(consensus)

        # Add domains to complete set of domains
        self.domain_set.update(domains.keys())

    def set_debug(self,val):
        self.debug=val
