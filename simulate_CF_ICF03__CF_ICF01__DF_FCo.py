'''
Created on Apr 6, 2015

@author: yanying
'''
import compare_CF_ICF03__CF_ICF01__DF_FCo as compare

import itertools as it


pathString = '/home/yanying/haihe_MultiComputer/ICF/ICF_journal_Revision/2ndRoundRevision/2ndRound_twoUserSimulation_Python/notes_TwoUserCase/plotsSimulation'

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }

# run #1:
run01 = False 
if run01:
    g13, g14, g23, g24 = 1, -1, 1, 1
    gParam = [g13, g14, g23, g24]
    
    # Note power Ps need to be greater than 0.5  
    PsList = [100]
    
    # consider Ps=100 as the bench mark 
    P3P4List = [[1,1], [10,10], [100,100], [1000,1000] ]
    
    numPoints=21
    
    simulations = []
    for (Ps,P3P4s) in it.product(PsList, P3P4List):
        P3, P4 = P3P4s
        print(Ps,P3,P4 , '--', g13, g14, g23, g24)
        
        oneSimulation = compare.Data(Ps, P3, P4, numPoints,
                             g13, g14, g23, g24,
                             font)
        
        powerParam = [Ps, P3, P4]
        # gParam = [g13, g14, g23, g24]
        simulations.append({'powerParam': powerParam, 'gParam': gParam, 'oneSimulation':oneSimulation})
        compare.display(oneSimulation, pathString, separate=False ,saveOnly=True)
        

# run #2:
run02 = False 
if run02:
    g13, g14, g23, g24 = 1, -1, 1, 1
    gParam = [g13, g14, g23, g24]
    
    # Note power Ps need to be greater than 0.5  
    PsList = [10]
    
    # consider Ps=100 as the bench mark 
    P3P4List = [[1,1], [10,10], [100,100], [1000,1000] ]
    
    numPoints=21
    
    simulations = []
    for (Ps,P3P4s) in it.product(PsList, P3P4List):
        P3, P4 = P3P4s
        print(Ps,P3,P4 , '--', g13, g14, g23, g24)
        
        oneSimulation = compare.Data(Ps, P3, P4, numPoints,
                             g13, g14, g23, g24,
                             font)
        
        powerParam = [Ps, P3, P4]
        # gParam = [g13, g14, g23, g24]
        simulations.append({'powerParam': powerParam, 'gParam': gParam, 'oneSimulation':oneSimulation})
        compare.display(oneSimulation, pathString, separate=False ,saveOnly=True)    
        
        
# run #3:
run03 = False 
if run03:
    g13, g14, g23, g24 = 2, -2, 2, 2
    gParam = [g13, g14, g23, g24]
    
    # Note power Ps need to be greater than 0.5  
    PsList = [10]
    
    # consider Ps=100 as the bench mark 
    P3P4List = [[1,1], [10,10], [100,100], [1000,1000] ]
    
    numPoints=21
    
    simulations = []
    for (Ps,P3P4s) in it.product(PsList, P3P4List):
        P3, P4 = P3P4s
        print(Ps,P3,P4 , '--', g13, g14, g23, g24)
        
        oneSimulation = compare.Data(Ps, P3, P4, numPoints,
                             g13, g14, g23, g24,
                             font)
        
        powerParam = [Ps, P3, P4]
        # gParam = [g13, g14, g23, g24]
        simulations.append({'powerParam': powerParam, 'gParam': gParam, 'oneSimulation':oneSimulation})
        compare.display(oneSimulation, pathString, separate=False ,saveOnly=True)                

# run #4:
# adjusted powers at the relays based on run #1
run04 = True 
if run04:
    g13, g14, g23, g24 = 1, -1, 1, 1
    gParam = [g13, g14, g23, g24]
    
    # Note power Ps need to be greater than 0.5  
    PsList = [100]
    
    # P3P4List = [[50,50], [60,60], [70,70], [80,80], [90,90] ]
    P3P4List = [[2000, 2000] ]
    
    # in my computation when P3>=2500, the first hop rate region will be completely
    # inside the second hop rate region.
    # P3P4List = [[2400,2400], [2500,2500], [2900,2900], [3900,3900], [4000,4000], [4900,4900] ]
    
    numPoints=21
    
    simulations = []
    for (Ps,P3P4s) in it.product(PsList, P3P4List):
        P3, P4 = P3P4s
        print(Ps,P3,P4 , '--', g13, g14, g23, g24)
        
        oneSimulation = compare.Data(Ps, P3, P4, numPoints,
                             g13, g14, g23, g24,
                             font)
        
        powerParam = [Ps, P3, P4]
        # gParam = [g13, g14, g23, g24]
        simulations.append({'powerParam': powerParam, 'gParam': gParam, 'oneSimulation':oneSimulation})
        compare.display(oneSimulation, pathString, separate=False ,saveOnly=True)
