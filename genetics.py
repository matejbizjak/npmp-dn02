from deap import creator, base, tools, algorithms
from numpy import random
from negative_core import get_model_params, model, simulate
from utils import eval_signal, oscilating, find_peaks

# CXPB  is the probability with which two individuals are crossed
# MUTPB is the probability for mutating an individual
CXPB = 0.5
MUTPB = 0.2
POPULATION_SIZE = 300
MAX_GENERATIONS = 15

def get_params():
    new_params = list(range(6))
    for i in range(0, 3):
        for j in range(0, 3):
            new_params[i] = random.randint(1,5)
            new_params[3 + j] = random.uniform(0.1, 1000)
    return new_params


def eval_one_max(individual):
    params = get_model_params()
    A, B, C = simulate(model, params[:4] + tuple(individual))

    peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = find_peaks(A)
    amplituda = 0
    perioda = 0

    if oscilating(peaks_vals_A, 0.01):
        amplituda, perioda = eval_signal(peaks_A, minimums_vals_A, peaks_vals_A)
        # print(amplituda, perioda)

    return [amplituda - perioda]

def genetic_algorithm():
    creator.create("FitnessMax", base.Fitness, weights=[1.0])
    creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()
    toolbox.register("individual", tools.initIterate, creator.Individual, get_params)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", eval_one_max)
    toolbox.register("mate", tools.cxTwoPoint)
    # toolbox.register("mate", tools.cxOnePoint)
    # toolbox.register("mate", tools.cxBlend, alpha=0.2)
    # toolbox.register("mutate", tools.mutFlipBit, indpb=0.2)
    # toolbox.register("mutate", tools.mutShuffleIndexes, indpb=0.1)
    # toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=0.5, indpb=0.5)
    toolbox.register("mutate", tools.mutUniformInt, indpb=0.5)
    toolbox.register("select", tools.selTournament, tournsize=3)

    random.seed(64)
    pop = toolbox.population(n=POPULATION_SIZE)

    print("Start of evolution")

    # Evaluate the entire population
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    print("  Evaluated %i individuals" % len(pop))

    # Extracting all the fitnesses of
    fits = [ind.fitness.values[0] for ind in pop]
    # Variable keeping track of the number of generations
    g = 0

    # Begin the evolution
    while max(fits) < 50 and g < MAX_GENERATIONS:
        # A new generation
        g = g + 1
        print("-- Generation %i --" % g)

        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))

        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            # cross two individuals with probability CXPB
            if random.random() < CXPB:
                toolbox.mate(child1, child2)
                # fitness values of the children
                # must be recalculated later
                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:
            # mutate an individual with probability MUTPB
            if random.random() < MUTPB:
                toolbox.mutate(mutant, [1,1,1,0,0,0],[4,4,4,1000,1000,1000])
                del mutant.fitness.values

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        print("  Evaluated %i individuals" % len(invalid_ind))
        # The population is entirely replaced by the offspring
        pop[:] = offspring
        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]

        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x * x for x in fits)
        std = abs(sum2 / length - mean ** 2) ** 0.5

        print("  Min %s" % min(fits))
        print("  Max %s" % max(fits))
        print("  Avg %s" % mean)
        print("  Std %s" % std)

        best_ind = tools.selBest(pop, 1)[0]
        print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))

    print("-- End of (successful) evolution --")

    best = tools.selBest(pop, 20)
    for k in range(10):
        print("%s. best individuals %s, %s" % (k, best[k], best[k].fitness.values))

    return [best[0], best[19]]
