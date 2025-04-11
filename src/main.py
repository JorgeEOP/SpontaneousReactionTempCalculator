from TemperatureCalculator import TemperatureCalculator as tc
import collections

TEMPERATURE_RANGE = range(1, 400)

if __name__ == '__main__':
    products: dict = collections.defaultdict(dict)
    reactants: dict = collections.defaultdict(dict)

    products["Fe2O3"] = 2
    reactants["Fe"] = 4
    reactants["O2"] = 3

    calculator = tc()
    calculator.setProducts(products)
    calculator.setReactants(reactants)
    calculator.setTemperatureRange(TEMPERATURE_RANGE)

    calculator.getGibbsEnergyAndTemperatures()