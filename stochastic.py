from __future__ import division
from collections import Counter
import random
import matplotlib.pyplot as plt
import numpy as np

##### PARAMETERS #####
'''
tissue_shape:           length of the square tissue patch, each subdivided unit
                        square can fit 1 bacterium
init_percent_full:      average percent of tissue that will have a bacterium at
                        initialization
subpop_types_count:     number of subpopulations included in the model. For
                        this model we will use 5 subpopulations as follows
                        * 0 : empty
                        * 1 : faecalibacterium prausnitzii
                        * 2 : bacteroides fragilis
                        * 3 : other commensal bacteria
                        * 4 : clostridium difficile (pathogenic bacteria)
init_subpop_cummprob:   probability of initializing into each of the subpops
growth_subpop_prob:     probability of each subpopulation dividing
death_subpop_prob:      probability of each subpopulation dying
antibiotic_eff:         efficacy of wide-range antibiotics (little prefernce
                        for bacterial domain, targets both commensal and
                        pathogenic bacteria)
immune_eff:             efficacy for hypersensitive immune system (overactive
                        immune system that targets commensal flora)
time_steps              number of timesteps to run
'''
tissue_shape = 5
subpop_types_count = 5
init_percent_full = 0.5
init_subpop_cummprob = [0, .1, .2, 1, 0] #probabilities of 0, .1, .1, .8, 0
growth_subpop_prob = [0] + ([.05] * (subpop_types_count - 1))
death_subpop_prob = [0] + ([.02] * (subpop_types_count - 1))
antibiotic_eff = 0.4
immune_eff = 0.2
time_steps = 1000

##### INITIALIZE #####
tissue_patch = [[]]

##### FUNCTIONALITY #####
def print_tissue():
    ''' Print the tissue patch '''
    global tissue_patch
    counter = Counter() #count bacteria in various subpopulations
    for subpop_num in range(subpop_types_count):
        counter[subpop_num] = 0
    for row_num in range(tissue_shape): #loop over tissue
        for col_num in range(tissue_shape):
            counter[tissue_patch[row_num][col_num]] += 1
    print "Subpop breakdown:\n|",
    for i in range(subpop_types_count):
        print "%i: %i |" % (i, counter[i]),
    print "\n\nTissue patch:"
    for row in tissue_patch: #print patch
        print row
    print

def heatmap():
    ''' Plot a heatmap of the tissue patch '''
    global tissue_patch
    plt.matshow(tissue_patch)
    plt.colorbar(ticks=range(subpop_types_count))
    plt.show()

def init_bacteria():
    ''' Initialize gut patch with some bacteria '''
    global tissue_patch
    tissue_patch = [[0 for i in range(tissue_shape)] for j in range(tissue_shape)]
    for row_num in range(tissue_shape): #loop over tissue
        for col_num in range(tissue_shape):
            if random.random() < init_percent_full:
                r = random.random()
                for k in range(tissue_shape):
                    if r < init_subpop_cummprob[k]:
                        tissue_patch[row_num][col_num] = k
                        break

def divide():
    ''' Bacteria can divide if there is an available spot '''
    global tissue_patch
    possible_transitions = [] #list possible neighbors
    for t1 in [-1,0,1]:
        for t2 in [-1,0,1]:
            if t1 != 0 or t2 != 0:
                possible_transitions.append((t1,t2))
    for row_num in range(tissue_shape): #loop over tissue
        for col_num in range(tissue_shape):
            curr = tissue_patch[row_num][col_num]
            for t1,t2 in possible_transitions:
                if random.random() < growth_subpop_prob[curr]:
                    try:
                         if  tissue_patch[row_num + t1][col_num + t2] == 0:
                             tissue_patch[row_num + t1][col_num + t2] = curr
                    except IndexError:
                        pass

def death():
    ''' Bacteria have a chance of dying at each timestep '''
    global tissue_patch
    for row_num in range(tissue_shape): #loop over tissue
        for col_num in range(tissue_shape):
            curr = tissue_patch[row_num][col_num]
            if curr != 0 and random.random() < death_subpop_prob[curr]:
                tissue_patch[row_num][col_num] = 0

def infect():
    ''' Infect the tissue patch with the C. Diff pathogen '''
    global tissue_patch
    m = tissue_shape // 2
    tissue_patch[m][m] = 4

def attack(antibiotic = True):
    ''' Attack bacteria in GI tract with wide-range antiobiotic (antibiotic =
    true) or with immune system (antibiotic = false)'''
    global tissue_patch
    for row_num in range(tissue_shape): #loop over tissue
        for col_num in range(tissue_shape):
            eff = antibiotic_eff if antibiotic else immune_eff
            if random.random() < eff:
                tissue_patch[row_num][col_num] = 0

def fecal_transplant():
    ''' Fecal bacteriotherapy transplant, transplant healthy bacteria from
    donor, supplementing commensal flora populations'''
    global tissue_patch
    for row_num in range(tissue_shape): #loop over tissue
        for col_num in range(tissue_shape):
            if random.random() < init_percent_full:
                r = random.random()
                for k in range(tissue_shape):
                    if r < init_subpop_cummprob[k]:
                        tissue_patch[row_num][col_num] = k
                        break

##### MAIN #####
if __name__ == "__main__":
    init_bacteria()
    print_tissue()
    divide()
    print_tissue()
    death()
    print_tissue()
    infect()
    print_tissue()
