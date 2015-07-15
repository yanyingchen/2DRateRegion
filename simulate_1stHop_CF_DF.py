'''
Created on Apr 6, 2015

@author: yanying
'''
import compare_1stHop_CF_DF as compare

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
    powers = [0.6, 1, 10, 100]
    
    numPoints=21
    
    simulations = []
    for Ps in powers:
        print(Ps , '--', g13, g14, g23, g24)
        
        oneSimulation = compare.Data(Ps, 
                             g13, g14, g23, g24,
                             font)
        simulations.append({'powerParam': Ps, 'gParam': gParam, 'oneSimulation':oneSimulation})
        compare.display(oneSimulation, pathString, saveOnly=True) 
    
    
# run #2: 
run02 = False
if run02:
    g13, g14, g23, g24 = 10, -10, 10, 10
    gParam = [g13, g14, g23, g24]
    
    # Note power Ps need to be greater than ( 1- 1/(100+100) )  
    powers = [0.996, 1, 10, 100]
    
    numPoints=21
    
    simulations = []
    for Ps in powers:
        print(Ps , '--', g13, g14, g23, g24)
        
        oneSimulation = compare.Data(Ps, 
                             g13, g14, g23, g24,
                             font)
        simulations.append({'powerParam': Ps, 'gParam': gParam, 'oneSimulation':oneSimulation})
        compare.display(oneSimulation, pathString, saveOnly=True)   
        
        
# run #3
run03 = False
if run03:
    g13, g14, g23, g24 = 2, -2, 2, 2
    gParam = [g13, g14, g23, g24]
    
    # Note power Ps need to be greater than ( 1- 1/(4+4) )  
    powers = [0.99, 1, 10, 100]
    
    numPoints=21
    
    simulations = []
    for Ps in powers:
        print(Ps , '--', g13, g14, g23, g24)
        
        oneSimulation = compare.Data(Ps, 
                             g13, g14, g23, g24,
                             font)
        simulations.append({'powerParam': Ps, 'gParam': gParam, 'oneSimulation':oneSimulation})
        compare.display(oneSimulation, pathString, saveOnly=True)  
     
# run #4     
run04 = False
if run04:
    gs =[[1,-1,1,1],
         [2, -2, 2, 2],
         [5, -5, 5, 5],
         [10, -10, 10, 10] ]
    
    # Note power Ps need to be greater than 
    #( 1- 1/(4+4) )  
    powers = [100,  #( 1- 1/(1+1) )
              25, #( 1- 1/(4+4) )    
              4, #( 1- 1/(25+25) )
              1 ] #( 1- 1/(100+100) )  
    
    
    numPoints=21
    
    simulations = []
    for i in range(len(powers)):
        Ps = powers[i]
        gParam = gs[i] 
        g13, g14, g23, g24 = gParam
        print(Ps, '--', g13, g14, g23, g24)
        
        oneSimulation = compare.Data(Ps, 
                             g13, g14, g23, g24,
                             font)
        simulations.append({'powerParam': Ps, 'gParam': gParam, 'oneSimulation':oneSimulation})
        compare.display(oneSimulation, pathString, saveOnly=True)      
    