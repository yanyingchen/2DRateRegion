'''
Created on May 5, 2015

@author: yanying
'''
import compareThreeICFSchemes as compare 


pathString = '/home/yanying/haihe_MultiComputer/research/ICF/journalRevision/2ndRoundRevision/simulation0405/finalPlots'

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
    
    P3P4s = [[20,20], [4,36] ]
    
    numPoints=21
    
    simulations = []
    for P3,P4 in P3P4s:
        print(P3,P4 , '--', g13, g14, g23, g24)
        
        oneSimulation = compare.Data(P3, P4, 
                                     font)
        simulations.append({'powerParam': [P3, P4], 'gParam': gParam, 'oneSimulation':oneSimulation})
        compare.display(oneSimulation, pathString, saveOnly=True) 
