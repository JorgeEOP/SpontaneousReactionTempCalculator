import collections
import math
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

# Units: g/mol
MoleculesMassDict = {
    "Fe": 55.84,
    "O2": 31.99,
    "Fe2O3": 159.68
}

# Units: J/K*g
MoleculesSpecificHeatCapacityDict = {
    "Fe": 12.40,
    "O2": 0, # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1
    "Fe2O3": -825.50
}



# TODO: to extend to more T ranges
class EntCalculatorBase:
    MolShomateCoeffDict = {
    "Fe": {
        "A": 23.97,
        "B": 8.36,
        "C": 0.00,
        "D": -0.00,
        "E": -0.00,
        "F": 0.26,
        "G": 62.06,
        "H": 7.78
    },
    "O2": {
        "A": 31.32,
        "B": -20.23,
        "C": 57.86,
        "D": -36.50,
        "E": -0.00,
        "F": -8.90,
        "G": 246.79,
        "H": 0.00
    },
    "Fe2O3": {
        "A": 93.43,
        "B": 108.35,
        "C": -50.86,
        "D": 25.58,
        "E": -1.61,
        "F": -863.20,
        "G": 161.07,
        "H": -825.50
        }
    }
    MoleculesStdFormationEnthalpy = {
        "Fe": 12.40,
        "O2": 0, # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1
        "Fe2O3": -822.2
    }

    def __init__(self):
        pass

    def getEnthalpy(self, molecule: str, temperatures: range):
        pass
    def getEntropy(self, molecule: str, temperatures: range):
        pass
