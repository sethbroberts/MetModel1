#script purpose: demo of using metmodel to find lethal reaction deletions

#import metmodel library
import metmodelCLI

#make a new model object
m = metmodelCLI.cb()	

#read in model from file
m.build_from_mm2('model_organisms/ssamodel2.txt', readquiet=True)	

#set objective
m.set_objective('Maximize', 'R_BIOMASS')	


#this variable will be either 'lethal' or 'nonlethal', depending on objective value
call = 'nonlethal'	

#solve the 'wild type model' (no deletions); verbose=False to supress output
m.solve(verbose=False)	

#print summary of wild type findings
print 'wild type' + '\t' + str(m.OBJECTIVE_VALUE) + '\t' + call		

for r in m.REACTIONS:
	#don't bother testing exchanges, biomass rxns
	if 'R_SRC' in r or 'R_ESC' in r or 'R_EXCH' in r or r == 'R_BIOMASS':		
		continue
	#lookup some information for this reaction
	name, reversible, notes, rawequation = m.REACTIONS[r]	
	#get the gpr for the reaction
	gpr = m.get_notes(r, 'Gene_association: ')				
	#get the reaction's subsystem (pathway(s))
	subsystem = m.get_notes(r, 'SUBSYSTEM: ')		
	#get the reaction equation		
	equation = m.get_equation(r)							

	#determine what 'default' upper and lower bounds are for each rxn
	if reversible:	
		#this is so you can reset these to default values after 'deleting' this rxn						
		default_lowerbound = -1000			
	else:
		default_lowerbound = 0
	default_upperbound = 1000
	
	#now delete this rxn by constraining it to zero flux
	m.set_constraint(r, 0, 0)				
	#solve the model
	m.solve(verbose=False)		
	#determine whether it's lethal or nonlethal by the objective value (biomass flux)			
	if float(m.OBJECTIVE_VALUE) < 1e-10:	
		call = 'lethal'
	else:
		call = 'nonlethal'
		
	#print results for this reaction	
	print r + '\t' + str(m.OBJECTIVE_VALUE) + '\t' + call + '\t' + gpr + '\t' + subsystem + '\t' + equation	
	
	#reset default upper and lower bounds, move to next rxn
	m.set_constraint(r, default_lowerbound, default_upperbound)		
		