import math
import collections
from EntCalculatorBase import EntCalculatorBase
import numpy as np

"""
Will get the Standard Entropy of a molecule:
  S° = A*ln(t) + B*t + C*t2/2 + D*t3/3 - E/(2*t2) + G
  Units: J/(mol*K)
"""
class EntropyCalculator(EntCalculatorBase):
    def __init__(self):
        super().__init__()
        pass
    
    def getEntropy(self, molecule: str, temperatures: range) -> dict:
        coefficients = self.MolShomateCoeffDict[molecule]
        stdEntropy =  collections.defaultdict(dict)

        print("Calculating Entropy of: {} in the range [{}, {})".format(molecule, temperatures.start, temperatures.stop))

        for t in temperatures:
            try:
                temp = t/1000
                entropy = coefficients["A"]* np.log(temp) + coefficients["B"]*temp + coefficients["C"]*math.pow(temp,2)/2 + \
                           coefficients["D"]*math.pow(temp,3)/3 - coefficients["E"]*(1/(2*math.pow(temp, 2))) + coefficients["G"]
                stdEntropy[t] = entropy

                #if (t > 294 and t < 300):
                #    print ("Temp: {} ; Standard entropy: {}".format(t, entropy) )

            except ZeroDivisionError:
                stdEntropy[t] = "NaN"
        print("---------------------------------------------")
        
        return stdEntropy