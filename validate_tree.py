''' 
Tests input sequence for whether it is a valid tree-sequence
'''

from parameter import ConsensusStatsCommandParser
from sequence import NeuriteSequence
import sys

version = 0.1

if __name__ == '__main__':
    try:
        sequence = sys.argv[1]

        a_stack = [[1,0,1]]
        for i in range(len(sequence)):
            char = sequence[i]
            #print(str(i)+" "+char)
            if len(a_stack) == 0:
                print("Sequence going beyond last complete tree at position "+str(i))
            if char == 'A':
                a_stack.append([0,1,0])
            elif char == 'C':
                a_stack[-1][a_stack[-1][0]+1] += 1
            elif char == 'T':
                a_stack[-1][a_stack[-1][0]+1] += 1
                while len(a_stack) > 0 and a_stack[-1][0] == 1:
                    top = a_stack.pop()
                    if top[1] > top[2]+1:
                        print("First tree larger than second tree ending at position "+str(i))
                    if len(a_stack) > 0:
                        a_stack[-1][a_stack[-1][0]+1] += top[1]+top[2]
                # After popping, remaining tree switches to second subtree
                if len(a_stack) > 0:
                    a_stack[-1][0] = 1

        if len(a_stack) > 0:
            print("Tree is incomplete")

    except (IOError, KeyboardInterrupt, IndexError) as e:
        print(str(e))
