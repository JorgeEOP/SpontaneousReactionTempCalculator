import collections
import math
from EntalpyCalculator import EnthalpyCalculator as hcalc

#Room Temp: 293.15K

"""
Dictionary with atom type and enthalpy
https://webbook.nist.gov/
"""

# Units: kJ/mol
MoleculesEnthalpySolidDict = {
    "Fe": 12.40,
    "O2": 0, # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1
    "Fe2O3": -825.50
}

# Units: J/K*mol
MoleculesEntropySolidDict = {
    "Fe": 34.76,
    "O2": 205.15,
    "Fe2O3": 87.28
}

d = collections.defaultdict(dict)

class TemperatureCalculator():
    """
    Class to calculate the temperature which will give a spontaneous Reaction
     Gibbs free energy: DeltaG = DeltaH - (T * DeltaS)
     For spontaneous reactions we need to have DeltaG < 0
    """
    __reactants: dict = collections.defaultdict(dict)
    __products: dict = collections.defaultdict(dict)
    __temperatureRange: range

    """
    Data model:
    reactants, products = {
            Molecule[String]: Ocurrence[integer],
    }
    """

    def __init__(self):
        self.__temperatureRange = range(250,300)
        pass

    def setReactants(self, reactants: dict) -> None:
        print("TemperatureCalculator::setReactants: ", reactants)
        self.__reactants = reactants

    def setProducts(self, products: dict) -> None:
        print("TemperatureCalculator::setProducts: ", products)
        self.__products = products
        
    def setTemperatureRange(self, _temperatureRange: range):
        self.__temperatureRange = _temperatureRange

    def getDeltaH(self) -> dict:
        """
        Calculates the change in enthalpy using the equation: deltaH = SumProducts(n * deltaH°_f) - SumReactants(m * deltaH°_f)
        Units: kJ/mol
        @return: Change in Enthalpy
        """
        deltaH: dict
        entalphyProducts: float = 0
        entalphyReactants: float = 0

        enthalpiesOfReactants = collections.defaultdict(dict)
        enthalpiesOfProducts = collections.defaultdict(dict)
        enthalpyCalc = hcalc()

        for reactant in self.__reactants:
            enthalpiesAndTemp = enthalpyCalc.getEnthalpy(reactant, range(0,700))
            print("Reactant:: ", reactant)
            for i in enthalpiesAndTemp:
                weightedEnthalpy = enthalpiesAndTemp[i] * self.__reactants[reactant]
                if (i > 295 and i < 300):
                    print(i, enthalpiesAndTemp[i], weightedEnthalpy)


        for product in self.__products:
            enthalpiesAndTemp = enthalpyCalc.getEnthalpy(product, range(0,700))
            print("product:: ", product)
            for i in enthalpiesAndTemp:
                weightedEnthalpy = enthalpiesAndTemp[i] * self.__products[product]
                if (i > 295 and i < 300):
                    print(i, enthalpiesAndTemp[i], weightedEnthalpy)
               
        """
        for product in self.__products:
            enthalpies = enthalpyCalc.getEnthalpy(product, range(0,700))
            enthalpiesOfProducts[product] = self.__products[product] * enthalpies
        """

        #print (enthalpiesOfReactants.keys())
        #print (enthalpiesOfProducts.keys())

        #for reactant in self.__reactants:
        #    print("Molecule: {} ; NumberOfMolecules: {} ; Enthalpy: {}".format(reactant, self.__reactants[reactant], MoleculesEnthalpySolidDict[reactant]))
        #    entalphyReactants += self.__reactants[reactant] * MoleculesEnthalpySolidDict[reactant]

        """
        for product in self.__products:
            print("Molecule: {} ; NumberOfMolecules: {} ; Enthalpy: {}".format(product, self.__products[product], MoleculesEnthalpySolidDict[product]))
            entalphyProducts += self.__products[product] * MoleculesEnthalpySolidDict[product]
            
        #deltaH = entalphyProducts - entalphyReactants
        print ("--------------------------------------------")

        return deltaH
        """

    def getDeltaS(self) -> float:
        """
        Calculates the change in entropy using the equation: deltaS = SumProducts(n * deltaS°_f) - SumReactants(m * deltaS°_f)
        Units: kJ/mol

        @return: Change in Entropy
        """
        deltaS: float
        entropyProducts: float = 0
        entropyReactants: float = 0

        for reactant in self.__reactants:
            print("Molecule: {} ; NumberOfMolecules {} ; Entropy: {}".format(reactant, self.__reactants[reactant], MoleculesEntropySolidDict[reactant]))
            entropyReactants += self.__reactants[reactant] * MoleculesEntropySolidDict[reactant]

        for product in self.__products:
            print("Molecule: {} ; NumberOfMolecules {} ; Entropy: {}".format(product, self.__products[product], MoleculesEntropySolidDict[product]))
            entropyProducts += self.__products[product] * MoleculesEntropySolidDict[product]



        deltaS = entropyProducts - entropyReactants
        print ("--------------------------------------------")

        return deltaS
    
    # T triggers a spontaneous Reaction when deltaG < 0
    def getGibbsEnergyAndTemperatures(self):
        try:
            deltaH: float = self.getDeltaH()
            deltaS: float = self.getDeltaS()
            deltaS = deltaS * (1000**-1)
            deltaG: float = 0
            print("deltaH: ", deltaH, "deltaS: ", deltaS)

            # Units: temperature has units of K
            for i, t in enumerate(TEMPERATURE_RANGE):
                deltaG = round(deltaH - t * deltaS, 2)
                print("Temperature: {} ; Gibbs Energy: {}".format(t, deltaG))


        except Exception as e:
            print("{}::Exception:: {}".format(self.getGibbsEnergyAndTemperatures.__name__, e))