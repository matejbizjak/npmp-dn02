from deap import creator, base, tools, algorithms
from numpy import random
from negative_core import get_model_params, model, simulate
from utils import eval_signal, oscilating, find_peaks

# CXPB  is the probability with which two individuals are crossed
# MUTPB is the probability for mutating an individual
# MUTFC is the mutation factor
CXPB = 0.5
MUTPB = 0.2
MUTFC = 0.5
POPULATION_SIZE = 300
MAX_GENERATIONS = 5


def get_six_params():
    new_params = list(range(6))
    for i in range(0, 3):
        for j in range(0, 3):
            new_params[i] = random.randint(1,5)
            new_params[3 + j] = random.uniform(0.1, 1000)
    return new_params

def get_ten_params():
    new_params = list(get_model_params())
    for i in range(0, 3):
        for j in range(0, 3):
            new_params[3 + i] = random.randint(1,5)
            new_params[6 + j] = random.uniform(0.1, 1000)
    return new_params


def eval_one_max(individual):
    subject = tuple(individual)
    if len(subject) == 6:
        params = get_model_params()
        A, B, C = simulate(model, params[:4] + subject)
    else:
        A, B, C = simulate(model, subject)

    peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = find_peaks(A)
    amplituda = 0
    perioda = 0

    if oscilating(peaks_vals_A, 0.01):
        amplituda, perioda = eval_signal(peaks_A, minimums_vals_A, peaks_vals_A)
        # print(amplituda, perioda)

    if perioda > 10:
        return[0]

    return [amplituda - perioda]

def mutate_ten_parameters(candidate):
    #(alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3)
    for idx, val in enumerate(candidate):
    	rnd = random.uniform(0, 1)
    	if rnd <= MUTPB:
            rnd2 = random.uniform(1 - MUTFC, 1 + MUTFC)
            new_value = val * rnd2
            if idx == 0:
                candidate[idx] = new_value % 10**5
                candidate[1] = 0.001 * candidate[idx] % 10**2
            elif idx == 1:
                candidate[idx] = new_value % 10**2
                candidate[0] = 1000 * candidate[idx] % 10**5
            elif idx == 2:
                candidate[idx] = new_value % 10**4
            elif idx == 3:
                candidate[idx] = int(new_value * rnd2)
            elif idx > 3 and idx < 7:
                candidate[idx] = int(new_value * rnd2) % 4
            else:
                candidate[idx] = new_value % 10 **3

    return candidate,

def genetic_algorithm(mode="six_params"):
    creator.create("FitnessMax", base.Fitness, weights=[1.0])
    creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()
    if mode == "all_params":
        toolbox.register("individual", tools.initIterate, creator.Individual, get_ten_params)
        toolbox.register("mutate", mutate_ten_parameters)
    else:
        toolbox.register("individual", tools.initIterate, creator.Individual, get_six_params)
        toolbox.register("mutate", tools.mutUniformInt, indpb=0.5)

    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", eval_one_max)
    toolbox.register("mate", tools.cxTwoPoint)
    # toolbox.register("mate", tools.cxOnePoint)
    # toolbox.register("mate", tools.cxBlend, alpha=0.2)
    # toolbox.register("mutate", tools.mutFlipBit, indpb=0.2)
    # toolbox.register("mutate", tools.mutShuffleIndexes, indpb=0.1)
    # toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=0.5, indpb=0.5)
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
    while max(fits) < 90 and g < MAX_GENERATIONS:
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
                if mode == "all_params":
                    toolbox.mutate(mutant)
                else:
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

    sorted = []
    for k in range(len(best[0])):
        if list(best[0])[k] > list(best[19])[k]:
            sorted.append([(best[19])[k], tuple(best[0])[k]])
        else:
            sorted.append([(best[0])[k], tuple(best[19])[k]])

    return sorted