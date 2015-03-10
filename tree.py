from sequence import NeuriteSequence

class TreeIndexLogic():
    ''' 
    A class which captures tree-specific logic given two characters.
    @param char1: First character
    @param char2: Second character
    '''
    def __init__(self, char1, char2):
        self.char1 = char1
        self.char2 = char2
        
    def get(self):
        if (self.char1 == 'A' and self.char2 == 'C') or (self.char2 == 'A' and self.char1 == 'C'):
            return 'A'
        if (self.char1 == 'A' and self.char2 == '-') or (self.char2 == 'A' and self.char1 == '-'):
            return 'A'
        if (self.char1 == 'T' and self.char2 == '-') or (self.char2 == 'T' and self.char1 == '-'):
            return 'T'
        if (self.char1 == 'C' and self.char2 == '-') or (self.char2 == 'C' and self.char1 == '-'):
            return 'C'
        raise Exception("Improper character alignment: "+self.char1+" with "+self.char2)

class TreeLogicFactory():
    '''
    Parses and processes the composite string to ultimately yield a single
    string which encapsulate the pairwise alignment.
    '''
    def __init__(self, str1, str2):
        self.str1 = str1
        self.str2 = str2
        
    def get_alignment(self):
        ''' 
        Simple function to merge two strings and produce a composite.
        @return: NeuriteSequence object representing the composite sequence.
        '''
        composite = ''
        for idx, char1 in enumerate(self.str1):
            char2 = self.str2[idx]
            if char1 == self.str2[idx]:
                composite += char1
            else:
                # Apply neuronal logic given two specific characters.
                composite += TreeIndexLogic(char1, char2).get()
        return composite
