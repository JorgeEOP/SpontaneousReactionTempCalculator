import collections
import math
from EntCalculatorBase import EntCalculatorBase
"""
https://de.webqc.org/molecular-weight-of-Fe.html

This class will calculate the standard Enthalpy of a molecule 
  H° = A*t + B*t2/2 + C*t3/3 + D*t4/4 - E/t + F - H

  Units : kJ/mol

  (((q = m * c * T)))

We need:
  mm: molecule Mass
  c: specific heat capacity
  t: temperature
"""
"""
MoleculesMassDict = {
    "Fe": 55.84,
    "O2": 31.99,
    "Fe2O3": 159.68
}
"""

class EnthalpyCalculator(EntCalculatorBase):
    def __init__(self):
        pass
    
    def getEnthalpy(self, molecule: str, temperatures: range) -> dict:
        coefficients = self.MolShomateCoeffDict[molecule]
        stdEntalpy =  collections.defaultdict(dict)

        print("Calculating Enthalpy of: {} in the range [{}, {})".format(molecule, temperatures.start, temperatures.stop))

        for t in temperatures:
            try:
                temp = t/1000
                enthalpy = coefficients["A"]*temp + \
                           coefficients["B"]*math.pow(temp,2)/2 +\
                           coefficients["C"]*math.pow(temp,3)/3 + \
                           coefficients["D"]*math.pow(temp,4)/4 - \
                           coefficients["E"]*(1/temp) + \
                           coefficients["F"] - \
                           coefficients["H"]
                stdEntalpy[t] = enthalpy + self.MoleculesStdFormationEnthalpy[molecule]

                #if (t > 294 and t < 300):
                if (t == 298):
                    print ("Temp: {} ; Standard enthalpy: {}".format(t, stdEntalpy[t]) )

            except ZeroDivisionError:
                stdEntalpy[t] = "NaN"

        print("---------------------------------------------")
        
        return stdEntalpy
    
        