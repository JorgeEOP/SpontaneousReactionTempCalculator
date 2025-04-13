import collections
import math
import pandas as pd
import numpy as np

from EnthalpyCalculator import EnthalpyCalculator as hcalc
from EntropyCalculator import EntropyCalculator as scalc

"""
Dictionary with atom type and enthalpy
https://webbook.nist.gov/
"""

class TemperatureCalculator():
    """
    Class to calculate the temperature which will give a spontaneous Reaction
     Gibbs free energy: DeltaG = DeltaH - (T * DeltaS)
     For spontaneous reactions we need to have DeltaG < 0
    """
    __reactants: dict = collections.defaultdict(dict)
    __products: dict = collections.defaultdict(dict)
    __temperatureRange: range = range(250,300)

    """
    Data model:
    reactants, products = {
            Molecule[String]: Ocurrence[integer],
    }
    """

    def __init__(self):
        pass

    def setReactants(self, reactants: dict) -> None:
        print("TemperatureCalculator::setReactants: ", reactants)
        self.__reactants = reactants

    def setProducts(self, products: dict) -> None:
        print("TemperatureCalculator::setProducts: ", products)
        self.__products = products
        
    def setTemperatureRange(self, _temperatureRange: range) -> None:
        self.__temperatureRange = _temperatureRange

    def getDeltaH(self) -> list[float]:
        """
        Calculates the change in enthalpy using the equation: deltaH = SumProducts(n * deltaH°_f) - SumReactants(m * deltaH°_f)
        Units: kJ/mol
        @return: Change in Enthalpy
        """
        deltaH: list[float]

        totalEnthalpyOfReactants = collections.defaultdict(dict)
        totalEnthalpyOfProducts = collections.defaultdict(dict)
        enthalpyCalc = hcalc()

        for reactant in self.__reactants:
            totalEnthalpyOfReactants[reactant] = enthalpyCalc.getEnthalpy(reactant, self.__temperatureRange)
        for product in self.__products:
            totalEnthalpyOfProducts[product] = enthalpyCalc.getEnthalpy(product, self.__temperatureRange)
                
        totalEnthalpyOfReactants = self.getWeightedEnergyReactant(totalEnthalpyOfReactants)
        totalEnthalpyOfProducts = self.getWeightedEnergyProduct(totalEnthalpyOfProducts)
        print("getDeltaH:: Weighted Enthalpy of 4Fe at 298K: {}".format(totalEnthalpyOfReactants["Fe"][298]))
        print("getDeltaH:: Weighted Enthalpy of 3O2 at 298K: {}".format(totalEnthalpyOfReactants["O2"][298]))
        print("getDeltaH:: Weighted Enthalpy of 2Fe2O3 at 298K: {}".format(totalEnthalpyOfProducts["Fe2O3"][298]))

        totalReactantsEnthalpy = self.getSumOfAllEnergies(totalEnthalpyOfReactants)
        totalProductsEnthalpy = self.getSumOfAllEnergies(totalEnthalpyOfProducts)
        print("getDeltaH:: Sum of Reactants Enthalpies at 298K: {}".format(totalReactantsEnthalpy[298]))
        print("getDeltaH:: Sum of Products Enthalpies at 298K: {}".format(totalProductsEnthalpy[298]))

        deltaH = [0.] * len(totalProductsEnthalpy)
        for i in range(len(totalProductsEnthalpy)):
            try:
                deltaH[i] = totalProductsEnthalpy[i] - totalReactantsEnthalpy[i]
            except:
                deltaH[i] = None
        print("getDeltaH:: Change in Enthalpy: {}".format(deltaH[298]))
        
        return deltaH
    
    def getDeltaS(self):
        """
        Calculates the change in entropy using the equation: deltaS = SumProducts(n * deltaS°_f) - SumReactants(m * deltaS°_f)
        Units: kJ/mol

        @return: Change in Entropy
        """
        deltaS: dict

        totalEntropyOfReactants = collections.defaultdict(dict)
        totalEntropyOfProducts = collections.defaultdict(dict)
        entropyCalc = scalc()

        for reactant in self.__reactants:
            totalEntropyOfReactants[reactant] = entropyCalc.getEntropy(reactant, self.__temperatureRange)
        for product in self.__products:
            totalEntropyOfProducts[product] = entropyCalc.getEntropy(product, self.__temperatureRange)
        
        totalEntropyOfReactants = self.getWeightedEnergyReactant(totalEntropyOfReactants)
        totalEntropyOfProducts = self.getWeightedEnergyProduct(totalEntropyOfProducts)

        totalEntropyOfReactants = self.getSumOfAllEnergies(totalEntropyOfReactants)
        totalEntropyOfProducts = self.getSumOfAllEnergies(totalEntropyOfProducts)

        deltaS = [0.] * len(totalEntropyOfProducts)
        for i in range(len(totalEntropyOfProducts)):
            try:
                deltaS[i] = totalEntropyOfProducts[i] - totalEntropyOfReactants[i]
            except:
                deltaS[i] = None

        return deltaS
    
    def getSumOfAllEnergies(self, allElementsAndEnergies: dict):
        sumOfAllEnergies = collections.defaultdict(dict)
        length = len(self.__temperatureRange) + 1
        result = [0.] * length

        # Loop over all molecules
        for molecule, tempAndEnergies in allElementsAndEnergies.items():
            for t, energy in tempAndEnergies.items():
                energies = tempAndEnergies.keys()
                try:
                    result[t] = result[t] + float(energy)
                except:
                    result[t] = energy

                
        return result
    
    def getWeightedEnergyReactant(self, allMoleculesEnergies):
        # Weight it by the number of molecules
        for mol, tempAndEnergy in allMoleculesEnergies.items():
            nMolecules = self.__reactants[mol]
            tempAndEnergy.update( (t, energy*nMolecules) for t, energy in tempAndEnergy.items() )

        return allMoleculesEnergies
    
    def getWeightedEnergyProduct(self, allMoleculesEnergies):
        # Weight it by the number of molecules
        for mol, tempAndEnergy in allMoleculesEnergies.items():
            nMolecules = self.__products[mol]
            #print("Total Energy of Product: {} with NMolexules: {}".format(mol, nMolecules))
            tempAndEnergy.update( (t, energy*nMolecules) for t, energy in tempAndEnergy.items() )
        return allMoleculesEnergies
    
    # T triggers a spontaneous Reaction when deltaG < 0
    def getGibbsEnergyAndTemperatures(self):
        try:
            deltaH = self.getDeltaH()
            print(2*"\n")
            deltaS = self.getDeltaS()
            print(2*"\n")

            deltaG = [0.] * len(deltaH)
            #print("deltaH: ", deltaH[298], "deltaS: ", deltaS[298])

            # Units: temperature has units of K
            for i, t in enumerate(self.__temperatureRange):
                dH = deltaH[i] * 1000 # kJ -> J
                dS = deltaS[i]
                tdS = t * dS
                deltaG[i] = round(dH - tdS, 2)
                if (t > 295 and t <= 300):
                    print("Temperature: {} \n" \
                          "dH: {} J/mol \n" \
                          "dH: {} kJ/mol \n" \
                          "dS: {} J/molK \n" \
                          "TdS: {} J/mol \n" \
                          "Gibbs Free Energy: {}".format(t, round(dH, 3), round(dH/1000, 3) , round(dS, 3), round(tdS, 2), round(deltaG[i], 2)))
                    print("---")
                #print("Temperature: {} ; dH: {} ; dS: {} ; TdS: {} ; Gibbs Free Energy: {}".format(t, dH, dS, tdS, deltaG[i]))
            
            #print("Temperature: {} ; Gibbs Energy: {}".format(298, deltaG[298]))

        except Exception as e:
            print("{}::Exception:: {}".format(self.getGibbsEnergyAndTemperatures.__name__, e))