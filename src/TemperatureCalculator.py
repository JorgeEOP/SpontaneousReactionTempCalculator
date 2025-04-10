import collections

#Room Temp: 293.15K

"""
Dictionary with atom type and enthalpy
all in kJ/mol
https://webbook.nist.gov/
"""

MoleculesEnthalpySolidDict = {
    "Fe": 12.40,
    "O2": 0, # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1
    "Fe2O3": -825.50
}

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

    def getDeltaH(self) -> float:
        """
        Calculates the change in enthalpy using the equation: deltaH = SumProducts(n * deltaH°_f) - SumReactants(m * deltaH°_f)
        @return: Change in Enthalpy
        """
        deltaH: float
        entalphyProducts: float = 0
        entalphyReactants: float = 0

        for product in self.__products:
            entalphyProducts += self.__products[product] * MoleculesEnthalpySolidDict[product]
        
        for reactant in self.__reactants:
            entalphyReactants += self.__reactants[reactant] * MoleculesEnthalpySolidDict[reactant]
            
        #print ("Products enthalpy: ", entalphyProducts)
        #print ("Reactant enthalpy: ", entalphyReactants)

        deltaH = entalphyProducts - entalphyReactants

        #Sum the elements in _reactants
        #Sum the elements in _products

        return deltaH

    def getDeltaS(self) -> float:
        """
        Calculates the change in enthalpy using the equation: deltaS = SumProducts(n * deltaS°_f) - SumReactants(m * deltaS°_f)
        @return: Change in Entropy
        """
        deltaS: float
        entropyProducts: float = 0
        entropyReactants: float = 0

        for product in self.__products:
            entropyProducts += self.__products[product] * MoleculesEntropySolidDict[product]
        for reactant in self.__reactants:
            entropyReactants += self.__reactants[reactant] * MoleculesEntropySolidDict[reactant]

        deltaS = entropyProducts - entropyReactants

        return deltaS
    
    