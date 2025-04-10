from TemperatureCalculator import TemperatureCalculator as tc
import collections

if __name__ == '__main__':
    products: dict = collections.defaultdict(dict)
    reactants: dict = collections.defaultdict(dict)

    products["Fe2O3"] = 2
    reactants["Fe"] = 4
    reactants["O2"] = 4

    calculator = tc()
    calculator.setProducts(reactants)
    calculator.setReactants(products)
    deltaH = calculator.getDeltaH()
    deltaS = calculator.getDeltaS()
    
    #print(deltaH)
    #print(deltaS)