import random
import operator
import math
import statistics
import matplotlib.pyplot as plt


def generate_chromosome():

    lower_limit = -5.12
    upper_limit = 5.12
    chromosome = []

    for i in range(10):
        gene = random.uniform(lower_limit, upper_limit)     
        gene = round(gene, 3)   
        chromosome.append(gene)
            
    return chromosome


def fitness_function(chromosome):
        
    summation = 0
    D = 10
        
    for i in range(10):
        
        square = pow(chromosome[i], 2)
        cos = math.cos( (2 * math.pi) * chromosome[i])
        cos_result = 10 * cos

        result = square - cos_result        
        summation += result

    fitness = (10*D) +  summation
    fitness = round(fitness, 3)
    return fitness


def fathers_selection(population_fitness):
    
    # Max value selection    
    max_value = max(population_fitness)    
    max_value = round(max_value)    
    
    # Adding percentage of maximum value
    PERCENTAGE = 25
    reference_value = (PERCENTAGE * max_value) / 100    
    reference_value = round(max_value + reference_value)    

    # Calculating slices based on fitness
    roulette_slices = []
    total_roulette_slices = 0
    for i in range(len(population_fitness)):
        
        slice_value = reference_value - population_fitness[i]
        slice_value = round(slice_value)
        total_roulette_slices += slice_value
        roulette_slices.append(slice_value)       

    # Selection of fathers with the roulette technique
    fathers = []    
    fathers_counter = 0

    # Execute according to the number of evaluations
    while fathers_counter < 2:
        
        # Generate a number between 0 and the sum of fitnesses of the population
        random_number = random.uniform(0, total_roulette_slices)

        individual_found = False
        individual_selected = 0
        i = 0
        slice_counter = roulette_slices[i]

        # Execute until we found a indivual position
        while individual_found == False:
            
            if slice_counter >= random_number:
                individual_selected = i
                individual_found = True
            else:
                i += 1
                slice_counter += roulette_slices[i]    
        
        if individual_selected not in fathers:
            fathers.append(individual_selected)
            fathers_counter += 1    
            
    return fathers
        

def fathers_mix(fathers, population):

    # Fathers assignment       
    first_father = population[fathers[0]]    
    second_father = population[fathers[1]]

    child = []
    for i in range(10):
        
        c_max = 0
        c_min = 0

        # Found upper limit and lower limit
        if first_father[i] >= second_father[i]:
            c_max = first_father[i]
            c_min = second_father[i]
        else:
            c_max = second_father[i]
            c_min = first_father[i]               

        # Calculating interval
        interval_value = c_max - (c_min)        

        ALPHA = .25
        exploration_factor = ALPHA * interval_value        

        # Calculating exploration interval
        exploration_upper = c_max + exploration_factor 
        exploration_lower = c_min - exploration_factor       

        # Correct gene selection (within limits)
        correct_gene = False
        while correct_gene == False:

            child_gene = random.uniform(exploration_lower, exploration_upper)
            child_gene = round(child_gene, 3)  

            if (child_gene >= -5.12) and (child_gene <= 5.12):            
                child.append(child_gene)
                correct_gene = True
    
    return child


def mutate_child(child, mutation_probability):
    
    new_mutate_child = []

    # Calculating the standard deviation
    standar_deviation = statistics.pstdev(child)

    for i in range(10):

        # Calculating the probability of mutation of the gene
        gene_probability = random.randint(1, 10)    
        is_mutated = False
        if gene_probability <= mutation_probability:
            is_mutated = True
        
        if is_mutated:        
            
            # Execute the process until a gene is found within the bounds
            correct_gene = False
            while correct_gene == False:

                # Calculating Gaussian distribution with a mean of 0
                gaussian_distribution = random.gauss(0, standar_deviation)

                # Calculating the new gene with the Gaussian distribution
                new_gene = child[i] + gaussian_distribution
                new_gene = round(new_gene, 3)  

                # Checking gene boundaries
                if (new_gene >= -5.12) and (new_gene <= 5.12):            
                    new_mutate_child.append(new_gene)
                    correct_gene = True
                    
        else:

            new_mutate_child.append(child[i])
    
    return new_mutate_child


def run():
    
    population_size = int(input('Tamaño de población: '))
    evaluations_number = int(input('Número de evaluaciones: '))
    mutation_probability = int(input('Probabilidad de mutación: ')) 

    evaluations_graph = []
    best_individual_graph = []    
    fitness_zero = False

    # Generate population
    population = []
    for i in range(population_size):
        population.append(generate_chromosome())
    
    # Execute according to the number of evaluations
    for i in range(evaluations_number):

        print("Evaluación: "+str(i))
        evaluations_graph.append(i)

        # Generate population fitness
        population_fitness = []
        for i in range(population_size):
            fitness = fitness_function(population[i])   
            population_fitness.append(fitness)                             

        # Generate decendents list minus one
        decendents_list = []
        for i in range(population_size-1):

            # Fathers selection
            fathers = fathers_selection(population_fitness)
            
            # Fathers mix
            child = []
            child = fathers_mix(fathers, population)      

            # Mutation        
            child_probability = random.randint(1, 10)    
            is_mutated = False
            if child_probability <= mutation_probability:
                is_mutated = True
            
            if is_mutated:                    
                new_child = mutate_child(child, mutation_probability)
                decendents_list.append(new_child)
            else:
                decendents_list.append(child)    

        # Order the population to found the best individual
        fitness_chromosome = {}

        for i in range(population_size):      
            aux = []
            aux.append(population_fitness[i])      
            aux.append(population[i])      
            fitness_chromosome[i] = aux
                    

        values_sort = sorted(fitness_chromosome.items(), key=operator.itemgetter(1), reverse=False)        

        best_individual = values_sort[0]
        best_individual = best_individual[1]
        best_individual = best_individual[1]                                            

        # Replace the population minus the best individual from the previous population
        population = []
        population.append(best_individual)
    
        for decendent in decendents_list:
            population.append(decendent)

        # Generate population fitness
        population_fitness = []
        for i in range(population_size):
            fitness = fitness_function(population[i])   
            population_fitness.append(fitness)     

        print("Mejor aptitud: "+ str(population_fitness[0]))  
        worst_individual = max(population_fitness)        
        print("Peor aptitud: "+str(worst_individual))        

        best_individual_graph.append(population_fitness[0])    
        
        # Evaluate if the algorithm should stop
        if population_fitness[0] < 1:
            fitness_zero = True
            break


    
    fig, ax = plt.subplots()
    ax.plot(evaluations_graph, best_individual_graph)
    plt.show()      

    if fitness_zero:
        print("Mejor aptitud menor a 1")


if __name__ == '__main__':
    run()
