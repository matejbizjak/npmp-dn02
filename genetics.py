from deap import creator, base, tools, algorithms
from numpy import random

from utils import eval_signal, oscilating, find_peaks
from pyDOE import lhs

# CXPB  is the probability with which two individuals are crossed
# MUTPB is the probability for mutating an individual
# MUTFC is the mutation factor
CXPB = 0.5
MUTPB = 0.2
MUTFC = 0.3
POPULATION_SIZE = 2000
MAX_GENERATIONS = 15


# model_type = 'positive'
model_type = 'negative'
_temp = __import__(model_type + '_core', globals(), locals(), ['get_model_params', 'model', 'simulate'], 0)
get_model_params = _temp.get_model_params
model = _temp.model
simulate = _temp.simulate

def get_ten_params():
    new_params = list(get_model_params())
    for i in range(0, 3):
        for j in range(0, 3):
            new_params[3 + i] = random.randint(1,5)
            new_params[6 + j] = random.uniform(0.1, 1000)
    return new_params

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
                new_value = (val + 1) % 4
                if new_value == 0:
                    new_value = 1
                candidate[idx] = new_value
            elif idx > 3 and idx < 7:
                new_value = (val + 1) % 4
                if new_value == 0:
                    new_value = 1
                candidate[idx] = new_value
            else:
                candidate[idx] = new_value % 10 **3

    return candidate,

def get_six_params():
    new_params = list(range(6))
    for i in range(0, 3):
        for j in range(0, 3):
            new_params[i] = random.randint(1,5)
            new_params[3 + j] = random.uniform(0.1, 1000)

    return new_params

def get_six_params_better():
    from app import denormalize
    m_range = [1, 4]
    K_range = [0.1, 1000]

    samples = lhs(6, samples=1, criterion='center')
    new_params = list(range(6))

    # denormalize
    new_params[0] = denormalize(samples[0][0], m_range[0], m_range[1], discrete=True)  # m1
    new_params[1] = denormalize(samples[0][1], m_range[0], m_range[1], discrete=True)  # m2
    new_params[2] = denormalize(samples[0][2], m_range[0], m_range[1], discrete=True)  # m3
    new_params[3] = denormalize(samples[0][3], K_range[0], K_range[1], discrete=False)  # K1
    new_params[4] = denormalize(samples[0][4], K_range[0], K_range[1], discrete=False)  # K2
    new_params[5] = denormalize(samples[0][5], K_range[0], K_range[1], discrete=False)  # K3
    return new_params

def eval_one_max(individual, params=get_model_params()):

    subject = tuple(individual)
    if len(subject) == 6:
        A, B, C = simulate(model, tuple(params[:4] + subject))
    else:
        A, B, C = simulate(model, subject)

    peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = find_peaks(A)

    if oscilating(peaks_vals_A, 0.01):
        amplituda, perioda = eval_signal(peaks_A, minimums_vals_A, peaks_vals_A)
        if amplituda > 0 and perioda > 0:
            return [amplituda / perioda]
    return[0]


def genetic_algorithm(mode="six_params"):
    model_params = get_model_params()
    creator.create("FitnessMax", base.Fitness, weights=[1.0])
    creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()
    if mode == "all_params":
        toolbox.register("individual", tools.initIterate, creator.Individual, get_ten_params)
        toolbox.register("mutate", mutate_ten_parameters)
    else:
        toolbox.register("individual", tools.initIterate, creator.Individual, get_six_params)
        # toolbox.register("individual", tools.initIterate, creator.Individual, get_six_params_better)
        toolbox.register("mutate", tools.mutUniformInt, indpb=0.2)

    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", eval_one_max, params=model_params)
    toolbox.register("mate", tools.cxTwoPoint)
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
    while g < MAX_GENERATIONS:
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
                    # toolbox.mutate(mutant, [1,1,1,0,0,0],[2,2,2,100,100,100])
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

    popul = list(pop)
    best = list(pop)
    print("len of pop %s:" % len(best))
    for k in range(len(popul)):
        if popul[k].fitness.values[0] < 10:
            best.remove(popul[k])

    best = [list(t) for t in set(tuple(element) for element in best)]
    print("len of best pop %s:" % len(best))

    for i in range(len(best)):
        A, B, C = simulate(model, get_model_params()[:4] + tuple(best[i]))

        peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = find_peaks(A)

        if oscilating(peaks_vals_A, 0.01):
            amplituda, perioda = eval_signal(peaks_A, minimums_vals_A, peaks_vals_A)
            print("a: ",amplituda,", p:",perioda)

    save_grapf(best)
    return best

def save_grapf(best):
    import numpy as nm
    import pandas as pd
    import seaborn as sb
    import matplotlib.pyplot as plt
    df = pd.DataFrame(nm.array(best), columns=["m1", "m2", "m3", "K1", "K2", "K3"])
    fig = sb.pairplot(df, markers="+")
    fig.savefig('plots/temp/genetics/' + 'population' + str(POPULATION_SIZE) +'_generations'
                + str(MAX_GENERATIONS) + '_samples' + str(len(best)) + model_type +'.png', format='png')
    plt.show()
