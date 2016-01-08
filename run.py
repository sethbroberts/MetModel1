#script purpose: demo of how to run metmodel

#load library
import metmodelCLI

#make a new constraint-based model instance
m = metmodelCLI.cb()

#build the model with info in text file(s)
m.build_from_mm2('model_organisms/cthmodel.txt')

#set the objective
m.set_objective('Maximize', 'R_BIOMASS')

#run FBA
m.solve()