'''
Created on Apr 6, 2015

@author: yanying
'''
import compare_2ndHop_ICF_FCo as compare 


pathString = '/home/yanying/haihe_MultiComputer/ICF/ICF_journal_Revision/2ndRoundRevision/2ndRound_twoUserSimulation_Python/notes_TwoUserCase/plotsSimulation'

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }

# run #1:
run01 = True 
if run01:
    g13, g14, g23, g24 = 1, -1, 1, 1
    gParam = [g13, g14, g23, g24]
    
    # consider Ps=100 as the bench mark 
    P3P4s = [[1,1], [10,10], [100,100], [1000,1000] ]
    
    numPoints=21
    
    simulations = []
    for P3,P4 in P3P4s:
        print(P3,P4 , '--', g13, g14, g23, g24)
        
        oneSimulation = compare.Data(P3, P4, 
                                     numPoints,
                                     font)
        simulations.append({'powerParam': [P3, P4], 'gParam': gParam, 'oneSimulation':oneSimulation})
        compare.display(oneSimulation, pathString, saveOnly=True) 
