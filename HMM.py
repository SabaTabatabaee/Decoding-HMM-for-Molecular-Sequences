from math import *
import sys


INF = float("inf")


class HMM:
    def __init__(self):
        self.alphabet = []
        self.eI = []      # background probability of each character in the alphabet
        self.eM = [{}]    # emission probability, one dictionary for each matching state in the model
        self.t = []       # transition probability, one dictionary for each set of states (D,M,I) in the model
        self.nstate = 0
    
    def load(self,hmmfile):
    # only load the first model in the given hmmfile if there are two or more models
        with open(hmmfile,'r') as fin:
            for line in fin:
                stream = line.strip().split()

                if stream[0] == "LENG":
                    self.nstate = int(stream[1])

                if stream[0] == "HMM": 
                    # read alphabet
                    self.alphabet = stream[1:]
                    
                    # read transition order
                    stream = fin.readline().strip().split()
                    trans_order = [(y[0]+y[3]).upper() for y in stream]
                    
                    # read the next line, if it is the COMPO line then ignore it and read one more line
                    stream = fin.readline().strip().split()
                    if stream[0] == "COMPO":
                        stream = fin.readline().strip().split()
                    
                    # now the stream should be at the I0 state; read the emission of the I0 state 
                    e = {}
                    for (x,y) in zip(self.alphabet,stream):
                        e[x] = -float(y)
                    self.eI.append(e)    
                    
                    # now the stream should be at the B state; read in the transition probability
                    stream = fin.readline().strip().split()
                    tB = {'MM':-INF,
                          'MD':-INF,
                          'MI':-INF,
                          'IM':-INF,
                          'II':-INF,
                          'ID':-INF,
                          'DM':-INF,
                          'DI':-INF,
                          'DD':-INF}
                    for x,y in zip(trans_order,stream):
                        tB[x] = -INF if y == '*' else -float(y)

                    self.t.append(tB)
                    break    

            for i in range(1, self.nstate+1):
                # read each set of three lines at a time
                stream = fin.readline().strip().split() # this one is the emission of the M state
                if float(stream[0]) != i:
                    print("Warning: incosistent state idexing in hmm file; expecting state " + str(i) + "; getting state " + stream[0])
                e = {}
                for x,y in zip(self.alphabet,stream[1:]):
                    e[x] = -INF if y == "*" else -float(y) 
                self.eM.append(e)    

                # the next line is the emission of the I state
                stream = fin.readline().strip().split() 
                e = {}
                for x,y in zip(self.alphabet,stream):
                    e[x] = -INF if y == "*" else -float(y) 
                self.eI.append(e)                    

                # the next line contains the transition probs
                stream = fin.readline().strip().split() # this one is the transition prob
                tB = {'MM':-INF,
                      'MD':-INF,
                      'MI':-INF,
                      'IM':-INF,
                      'II':-INF,
                      'ID':-INF,
                      'DM':-INF,
                      'DI':-INF,
                      'DD':-INF}
                for x,y in zip(trans_order,stream):
                    tB[x] = -INF if y == '*' else -float(y)
                self.t.append(tB)      

    def compute_llh(self, query):
        # Compute likelihood of an aligned query to the HMM
        # Return -inf if the query is not properly aligned
        j = 0
        prev_state = 'M'
        llh = 0

        for c in query:
            if c == '-': # gap -> deletion
                curr_state = 'D'
            elif c >= 'a' and c <= 'z': # lowercase -> insertion
                curr_state = 'I'
            elif c >= 'A' and c <= 'Z': # uppercase -> match
                curr_state = 'M'
            else: # encounter invalid symbol
                #print("encountered invalid symbol!")
                return -INF

            trans = prev_state + curr_state

            # penalize the transition
            if trans in self.t[j]:
                #print(self.t[j][trans])
                llh += self.t[j][trans]
            else: # encounter invalid transition
                #print("invalid transition  " + trans)
                return -INF 
            
            # transit: update the state index
            j += (curr_state != 'I') # move to the nex state unless this is moving towards an 'I'
            if j > self.nstate: # reach end of the HMM chain but not end of the sequence
                #print("reach end of the HMM chain but not end of the sequence")
                return -INF

            # penalize the emission
            if curr_state == 'M':
                #print(self.eM[j][c])
                llh += self.eM[j][c]
            elif curr_state == 'I': 
                #print(self.eI[j][c.upper()])   
                llh += self.eI[j][c.upper()]
            
            # update state
            prev_state = curr_state

        if j != self.nstate: # does not reach the E state at the end
            #print("does not reach the E state at the end")
            return float("-inf")
        trans = prev_state + 'M'
        llh += self.t[j][trans]    
        return llh

    def Viterbi(self, query):
        L = len(query) 
        K = self.nstate 
                # Initialize Viterbi matrices
        VM = [[-INF] * (L + 1) for _ in range(K + 1)]
        VI = [[-INF] * (L + 1) for _ in range(K + 1)]
        VD = [[-INF] * (L + 1) for _ in range(K + 1)]
        # Initialize traceback matrices
        traceback_M = [[-1] * (L + 1) for _ in range(K + 1)]
        traceback_I = [[-1] * (L + 1) for _ in range(K + 1)]
        traceback_D = [[-1] * (L + 1) for _ in range(K + 1)]

        VM[0][0] = 0  
                # Fill in the Viterbi matrices
        for i in range(1, L + 1):
            for j in range(K+1):
                if j >0:
                    match_emission = self.eM[j].get(query[i - 1], -INF)
                    max_VM, prev_state = max(
                        (VM[j - 1][i-1 ] + self.t[j-1]['MM'], 'M'),
                        (VI[j - 1][i-1] + self.t[j-1]['IM'], 'I'),
                        (VD[j - 1][i-1] + self.t[j-1]['DM'], 'D'),
                        key=lambda x: x[0]
                    )
                    VM[j][i] = match_emission + max_VM
                    traceback_M[j][i] = prev_state

                # Compute VI[j][i]
                match_emission = self.eI[j].get(query[i -1], -INF)
                max_VI, prev_state = max(
                    (VM[j][i - 1] + self.t[j]['MI'], 'M'),
                    (VI[j][i - 1] + self.t[j]['II'], 'I'),
                    (VD[j][i - 1] + self.t[j]['DI'], 'D'),
                    key=lambda x: x[0]
                )
                VI[j][i] = match_emission + max_VI
                traceback_I[j][i] = prev_state

                # Compute VD[j][i]
                if j > 0:
                    max_VD, prev_state = max(
                        (VM[j - 1][i] + self.t[j - 1]['MD'], 'M'),
                        (VI[j - 1][i] + self.t[j - 1]['ID'], 'I'),
                        (VD[j - 1][i] + self.t[j - 1]['DD'], 'D'),
                        key=lambda x: x[0]
                    )
                    VD[j][i] = max_VD
                    traceback_D[j][i] = prev_state


        # Backtrack to find Vscore, aln
        
        alignment = []
        if K<=L:
            j = K
            Vscore, best_state = max(
            (VM[K][L], 'M'),
            (VI[K][L], 'I'),
            (VD[K][L], 'D'),
            key=lambda x: x[0]
            )
            state = best_state
            for i in range(L, 0, -1):
                if state == 'M':
                    alignment.append(query[i - 1].upper())
                    state = traceback_M[j][i]
                    j -= 1
                elif state == 'I':
                    alignment.append(query[i - 1].lower())
                    state = traceback_I[j][i]
                    
                else:
                    alignment.append('-')
                    state = traceback_D[j][i]
                    
        if K>L:
            Vscore, best_state = (VD[K][L], 'D')
            state = best_state
            j=K
            i=L
            for _ in range(K):
                if state == 'M':
                    alignment.append(query[i - 1].upper())
                    state = traceback_M[j][i]
                    j=j-1
                    i=i-1
                    
                elif state == 'I':
                    alignment.append(query[i - 1].lower())
                    state = traceback_I[j][i]
                    j=j-1
                else:
                    alignment.append('-')
                    state = traceback_D[j][i]
                    j=j-1
                
                    
        alignment.reverse()   
        aln = ''.join(alignment)
        if Vscore == float("-inf"):
            aln = '$'
        return Vscore, aln

    




        

        


        






        
    
