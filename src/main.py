from TemperatureCalculator import TemperatureCalculator as tc
from EntropyCalculator import EntropyCalculator as scalc
import collections

TEMPERATURE_RANGE = range(0, 700)

if __name__ == '__main__':
    products: dict = collections.defaultdict(dict)
    reactants: dict = collections.defaultdict(dict)

    products["Fe2O3"] = 2
    reactants["Fe"] = 4
    reactants["O2"] = 3

    calculator = tc()
    calculator.setProducts(products)
    calculator.setReactants(reactants)
    #calculator.setTemperatureRange(TEMPERATURE_RANGE)

    calculator.getDeltaH()
    #calculator.getGibbsEnergyAndTemperatures()

    #entrpopyCalc = scalc()
    #entrpopyCalc.getEntropy("Fe", range(0,700))