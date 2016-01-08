#script purpose: make and analyze metabolic models in Python
#cb is class for constraint-based models 
	#uses Python 2.4.2 with standard libraries
	#works with Python 2.7.10, again with standard libraries
	#requires: glpsol (from glpk, https://www.gnu.org/software/glpk/)
	#also uses eq_current.py module, written to deal with parsing reaction equations, metabolites, compartments, etc.
	#this version omits mapGPR.py module, 
	#   written to read / parse / evaluate boolean GPR statements, etc.

"""
data structures for a constraint-based model:

modelid = modelid (string)

modelname = modelname (string)

compartments -> { compartmentname : {id:compartmentname, outside:outsidecompartment(None)}, ... }

species -> { speciesid : {id:speciesid, name:speciesname, charge:speciescharge, compartment:speciescompartment, boundarycondition:speciesBC}, ... }

reactions -> { reactionid: ( name, rev, {<notes>}, [[(r1, coef1), (r2, coef2), ... ], [(p1, coef1), (p2, coef2)]] ), ... }

NOTES:
1. currently, reversibility is determined by parsing rxnequation when reading tab-delimited input files,
	i.e., the column 'REVERSIBILITY' is ignored. Might eventually change this, perhaps eliminate column from input, or use as a check.
"""

import os, re, time, pickle		#standard Python modules
import eq_current				#custom Python module


#regular expression to capture ec numbers
ecnum_re = re.compile(r"""\d+\.(\d|-)+\.(\d|-)+\.(\d|-)+""")

#these are reactions where the equation is not quite the same across the palsson models for all organisms'
DISCREPANCIES = {'R_DHFS':1, 'R_DHPS2':1, 'R_MTHFR2':1, 'R_MTHFR3':1, 'R_METS':1, 
				'R_MTHFCm':1, 'R_GTPCI':1, 'R_MTHFD':1, 'R_MTHFC':1, 'R_MTHFD2':1, 
				'R_MTHFDm':1, 'R_QULNS':1}

#dictionary mapping one letter abbreviation used as suffix on species ID to corresponding compartment				
abbrev2compartment = {
						'c':('Cytosol', 'Extraorganism'),
						'r':('EndoplasmicReticulum', 'Cytosol'),
						'e':('Extraorganism', False),
						'g':('GolgiApparatus', 'Cytosol'),
						#'l':('Lysosome', 'Cytosol'),
						'l':('Reservosome', 'Cytosol'),
						'm':('Mitochondria', 'Cytosol'),
						'n':('Nucleus', 'Cytosol'),
						#'x':('Peroxisome', 'Cytosol'),
						#'y':('Glycosome', 'Cytosol')
						'x':('Glycosome', 'Cytosol'),
						'v':('Vacuole', 'Cytosol'),
						'b':('Extraorganism', False),
						'p':('BtwnMitoInnerOuter', 'Cytosol')
				
											}
											
def derive_coef (original_raw):
	#Given coeficient from eq data structure, derive the coeficient to be used in the *.lp file
	coef = ' '
	raw = str(original_raw)
	if raw[0] == '-':
		integ = raw[1:]
		if integ == '1':
			coef = coef + '-'
		else:
			coef = coef + '-' + integ
	else:
		integ = raw
		if integ == '1':
			coef = coef + '+'
		else:
			coef = coef + '+' + integ
	return coef
											
def ensure_boolean (id, val):
	#ensure that boundaryCondition is a boolean value
	if type(val) == type(True):
		return val
	else:
		assert (('rue' in val) or ('alse' in val)), 'Something wrong with boolean value: %s %s' % (id, val)
		if 'rue' in val:
			boolvar = True
		elif 'alse' in val:
			boolvar = False
		return boolvar
		

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class cb:
	"""
	Class for creation and manipulation of constraint-based metabolic models.
		
	"""
	
	def __init__ (self):
	
		#make a timestamp and use it as default 'model ID'...
		timestamp = time.strftime("%Y_%m_%d_%H_%M_%S")
		self.MODEL_ID = 'Model' + timestamp

		self.MODEL_NAME = ''
		self.COMPARTMENTS = {}
		self.SPECIES = {}
		self.REACTIONS = {}
		self.SOURCES = []
		self.ESCAPES = []
		self.EXCHANGES = []
		self.NOTSOURCES = []
		self.NOTESCAPES = []
		self.OBJECTIVE_EQUATION = ''
		self.GENES, self.TRANSCR, self.PROTS, self.REACTS, self.COMPLEXES, self.ISOZYMES =	{}, {}, {}, {}, {}, {}
		self.PROTEIN2GENE = {}
		self.SIMPLEGPR = {}

				
		#default max/min value for fluxes
		self.VMAX = '1000'
		
		#holds user-specified reaction flux constraints
		self.CONSTRAINTS = {}
		
		self.OBJECTIVE = ('Maximize', 'R_biomass')
		self.STATUS = ''
		self.OBJECTIVE_VALUE = ''
		self.REACTION2FLUXVALUE = {}
		self.MINBIOMASS = '0.001'			
				
								
	def set_id (self, ID):
		"Sets the name of the model. For now, model ID (as opposed to model name) will be a timestamp."
		self.MODEL_ID = ID


	def set_name (self, ID):
		"Sets the name of the model. For now, model ID (as opposed to model name) will be a timestamp."
		self.MODEL_NAME = ID
		

	def set_objective (self, goal, ID):
		"Set whether to maximize/minimize and which reaction. Example: m.set_objective('Maximize', 'R_biomass')."
		self.OBJECTIVE = (goal, ID)
	
	
	def delete_reaction (self, id):
		"Given a reactionID, delete this key, value pair from REACTIONS. Does not delete reaction species from SPECIES."
		if id in self.REACTIONS:
			del self.REACTIONS[id]
		else:
			print 'WARNING--cannot delete %s: not in REACTIONS' % (id)
			
	
	def get_equation (self, id):
		"Given a reactionID, prettyprint the reaction equation."
		name, reversible, notes, equation = self.REACTIONS[id]
		reactionequation = eq_current.makestring(equation, reversible)
		return reactionequation
		
		
	def add_species (self, id, name, compartment, charge, boundaryCondition):
		"Write a new species into the species list for the model. Example: m.add_species('M_h2o_c', 'water', 'cytosol', '0', 'False'). Params are ID, name, compartment, charge, boundaryCondition."
		self.SPECIES[id] = {'id':id, 'name':name, 'compartment':compartment, 'charge':charge, 'boundaryCondition':boundaryCondition}


	def add_compartment (self, id, outside=None):
		"Add a new compartment to the model. Example: m.add_compartment('Cytosol', outside='Extraorganism')"
		self.COMPARTMENTS[id] = {'id':id}
		if outside:
			self.COMPARTMENTS[id]['outside'] = outside
		
		
	def set_constraint (self, id, lbound, ubound):
		"Given a reaction ID, set lbound and ubound. Example: m.set_constraint('R_UNK2', 0, 1000)."
		if id in self.REACTIONS:
			self.CONSTRAINTS[id] = (str(lbound), str(ubound))
		else:
			print 'WARNING--cannot set constraint for %s: not in REACTIONS' % (id)
	
	
	def unset_constraint (self, id):
		"Given a reaction ID, reset lbound and ubound to defaults. Example: m.unset_constraint('R_UNK2')."
		if id in self.REACTIONS and id in self.CONSTRAINTS:
			del self.CONSTRAINTS[id]
		elif not id in self.REACTIONS:
			print 'WARNING--cannot unset constraint for %s: not in REACTIONS' % (id)
	
	
	def print_constraints (self):
		"Print all constraints."
		orderedc = self.CONSTRAINTS.keys()
		orderedc.sort()
		for id in orderedc:
			print id + '\t' + self.CONSTRAINTS[id][0] + '\t' + self.CONSTRAINTS[id][1]
	
	
	def reset_vmax (self, newvalue):
		"Resets the default flux limits to the provided value."
		newvalue_str = str(newvalue)
		for constraint in self.CONSTRAINTS:
			lbound, ubound = self.CONSTRAINTS[constraint]
			if lbound == '-' + self.VMAX:
				lbound = '-' + newvalue_str
			if ubound == self.VMAX:
				ubound = newvalue_str
			self.CONSTRAINTS[constraint] = (lbound, ubound)
		self.VMAX = newvalue_str
				
				
	def	add_note (self, ID, notetext):
		"Add a note to the corresponding reaction. Example: m.add_note('R_PGM', 'SUBSYSTEM: ss glycolysis')."
		if not ID in self.REACTIONS:
			pass
			#print 'WARNING--cannot add "%s" to notes for %s: %s not in REACTIONS' % (notetext, ID, ID)
		else:	
			name, rev, notes, eq = self.REACTIONS[ID]
			notes[notetext] = 1
			self.REACTIONS[ID] = (name, rev, notes, eq)


	def	delete_note (self, ID, notetext):
		"Delete a note from the corresponding reaction. Example: m.delete_note('R_PGM', 'SUBSYSTEM: ss glycolysis')."
		if not ID in self.REACTIONS:
			print 'WARNING--cannot delete "%s" from notes of %s: %s not in REACTIONS' % (notetext, ID, ID)
		else:	
			name, rev, notes, eq = self.REACTIONS[ID]
			if notetext in notes:
				del notes[notetext]
				self.REACTIONS[ID] = (name, rev, notes, eq)
			else:
				print 'WARNING--cannot delete "%s" from notes of %s: %s not in notes' % (notetext, ID, notetext)


	def get_notes (self, reaction, tagstring):
		"Get specific categories of information from notes; categories are indicated by tagstring, e.g., 'SUBSYSTEM: ' for pathway info."
		results = {}
		name, rev, notes, eq = self.REACTIONS[reaction]
		for note in notes:
			if tagstring in note:
				results[note[len(tagstring):]] = 1
		if results == {}:
			results = {'.':1}
		return (' ').join(results.keys())


	def add_reaction (self, ID, name, rev, notes, equation):
		"Add a new reaction into the model. Example: m.make_reaction('R_ss_biomass', 'ssa biomass', 'false', {'CONFIDENCE: 1':1, 'SUBSYSTEM: biomass':1, 'GPR: ':1, 'EC Number: ':1}, [[('M_atp_c', '1')],[('M_adp_c', '1'), ('M_pi_c', '1')]])"
		if ID in self.REACTIONS:
			print ID, 'already in REACTIONS'
		else:
			self.REACTIONS[ID] = (name, rev, notes, equation)
			if ID in DISCREPANCIES:
				warning_equation = eq_current.makestring(equation, rev)
				#print ID, 'discrepant across models. Using:', warning_equation
			for species, coef in (equation[0] + equation[1]):
				#add compartment if necessary
				compartment, outside = abbrev2compartment[species[-1:]]
				if not compartment in self.COMPARTMENTS:
					cb.add_compartment(self, compartment, outside)
				#add species if necessary
				if not species in self.SPECIES:
					name, db_compartment, charge, boundaryCondition = '.', '.', '.', 'false'
					cb.add_species(self, species, name, compartment, charge, boundaryCondition)
				

	def write_lp (self, lpfilename):
		"Write current model in *.lp file format. Provide a name for the file. Automatically called by 'solve' method."
		outfile = open(lpfilename, 'w')
		
		mets = {}
		constraints = {}

		for ID in self.REACTIONS:
			
			name, reversible, notes, equation = self.REACTIONS[ID]
						
			#if there have been specific constraints already set, use them...
			if ID in self.CONSTRAINTS:
				lbound, ubound = self.CONSTRAINTS[ID]
			#otherwise create constraints on the fly, using reaction reversibility and self.VMAX...
			else:
				ubound = self.VMAX
				if bool(reversible):
					lbound = '-' + self.VMAX
				else:
					lbound = '0'
			constraints[lbound + ' <= ' + ID + ' <= ' + ubound] = 1
				
			#create a data structure called 'mets': keys are metabolites, values are (reactionID, coef); this is essentially "S * v"
			for reactant in equation[0]:
				species, coef = reactant[0], '-' + str(reactant[1])
				if '_b' == species[-2:]:
					continue
				if not species in mets:
					mets[species] = []
				mets[species].append((ID, coef))
				
			for product in equation[1]:
				species, coef = product[0], product[1]
				if '_b' == species[-2:]:
					continue
				if not species in mets:
					mets[species] = []
				mets[species].append((ID, coef))

		#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		#WRITE *.LP FILE....

		assert self.OBJECTIVE[1], 'No objective has been defined.'
		if not self.OBJECTIVE[1] in self.REACTIONS:
			print
			print '!! Warning:', self.OBJECTIVE[1], 'is not a reaction in the model !!'
			print

		#print .lp file header
		print >>outfile, '\\\ '
		print >>outfile, '\n'
		print >>outfile, '\\\ Objective function'
		print >>outfile, self.OBJECTIVE[0]
		print >>outfile, '  Z : ' + self.OBJECTIVE[1]
		print >>outfile, '\n'
				
		#print out 'Subject To' part of .lp file
		print >>outfile, '\\\ Mass balance equations'
		print >>outfile, 'Subject To'		

		for m in mets:
			fcs = mets[m]
			line = '  ' + m + ' :'
			
			for fc in fcs:
				flux = fc[0]
				coef = derive_coef(fc[1])
				line = line + coef + ' ' + flux
			line = line + ' = 0'
			print >>outfile, line

		#print constraints
		print >>outfile, '\n'
		print >>outfile, '\\\ Flux constraints'
		print >>outfile, 'Bounds'
		print >>outfile, '\n'

		for c in constraints:
			print >>outfile, '  ' + c

		#print *.lp file suffix...
		print >>outfile, '\n'
		print >>outfile, 'End'
		
		
	def solve (self, out=False, verbose=True):
		"Run glpsol to see if solution exists. Argument is out=<fn> (if no filename given, just solves without writing output to a file, for checking purposes)."
		
		#if no escapes have been specified, make escapes on all metabolites in the model
		if self.ESCAPES == [] and self.EXCHANGES == []:
			print '# No escapes currently specified. Adding escape fluxes to all metabolites in model.'
			cb.set_escapes(self, self.SPECIES.keys())

		#make timestamp...
		timestamp = time.strftime("%Y_%m_%d_%H_%M_%S")
			
		if out:
			#set names of outputfiles
			lpfilename = out + '.' + timestamp + '.lp'
			rawoutfilename = out + '.' + timestamp + '.out'
			xlsfilename = out + '.' + timestamp + '.xls'
			
		else:
			#if out not specified, make tmp filenames (these files deleted below in this case)
			lpfilename = 'tmp.' + timestamp + '.lp'
			rawoutfilename = 'tmp.' + timestamp + '.out'
			xlsfilename = 'tmp.' + timestamp + '.xls'
			
		#write the *.lp file
		cb.write_lp(self, lpfilename)

		#construct glpsol command and execute, following calls glpsol from .lib; original command commented out below
		#command = '/Users/seth/.lib/python/glpsol --cpxlp ' + lpfilename + ' -o ' + rawoutfilename + ' > glpsol.log'
		command = 'glpsol --cpxlp ' + lpfilename + ' -o ' + rawoutfilename + ' > glpsol.log'
		os.system(command)
		
		self.REACTION2FLUXVALUE = {}
		
		#read rawoutput file, parse results
		file = open(rawoutfilename)
		lines = file.readlines()
		for i, line in enumerate(lines):
			if line == '': break
			line = line.rstrip()
			if line == '': continue
			col = line.split()
			
			#collects info on whether optimization was OK, and if so, what was the value of the Objective fnc
			if 'Status:' == line[0:7]:
				tmp = line.split()
				self.STATUS = tmp[1]
			if 'Objective:' == line[0:10]:
				tmp = line.split()
				self.OBJECTIVE_VALUE = float(tmp[3])
				
			#skipping irrelevant lines	
			if len(col) > 1:
				#find lines that have the name of a flux as col[1]
				if col[1] in self.REACTIONS:
					if len(col) < 4:
						nextline = lines[i+1]
						nextcol = nextline.split()
						self.REACTION2FLUXVALUE[col[1]] = nextcol[1]
					else:
						self.REACTION2FLUXVALUE[col[1]] = col[3]
							
		
		#send results to *.xls file
		if out:
			cb.list_reactions(self, out=xlsfilename, showfluxvalues=True)

		#if you didn't ask to save the output, delete the tmp output files
		if not out:
		
			if verbose:
				cb.list_reactions(self, showfluxvalues=True)
			
			command1 = 'rm ' + lpfilename
			command2 = 'rm ' + rawoutfilename

			os.system(command1)
			os.system(command2)
				
		
	def list_reactions (self, out=False, showfluxvalues=True):
		"Prints a list of reactions from current model, organized by path, then ecnumber. Arguments are out=<fn>, showfluxvalues=<True/False>. Defaults are False, True."
		cache = {}
		for reaction in self.REACTIONS:
			name, reversible, notes, equation = self.REACTIONS[reaction]
			reactionequation = eq_current.makestring(equation, reversible)
			#search for any ec numbers and pathways in reaction notes; it IS possible for there to be > 1 ec or pathway for a given reaction
			confidence, gpr = '?', '?'
			holder = {'pathways':{}, 'ecs':{}}
			ref, prr = '.', '.'
			for note in notes:
				if 'SUBSYSTEM: ' in note:
					holder['pathways'][note[11:]] = 1
				if 'EC: ' in note:
					holder['ecs'][note[4:]] = 1
				if 'CONFIDENCE: ' in note:
					confidence = note[12:]
				if 'GPR: ' in note:
					gpr = note[5:]
				if 'Protein_reaction_relation: ' in note:
					prr = note[note.find(' == ') + 4:]
				if 'PMID: ' in note: # and not 'review' in note and not 'related_organism' in note 
					ref = ref + note.split(',')[0][6:] + ' '
			if not ref == '.':
				ref = 'PMIDs: ' + ref[1:-1]
			
			#construct grr
			grr = prr
			for i in prr.split():
				if '(' in i:
					i = i[1:]
				if ')' in i:
					i = i[:-1]
				if i in self.PROTEIN2GENE:
					grr = grr.replace(i, self.PROTEIN2GENE[i])
					
			#add to cache for printing...																		
			for pathwayname in holder['pathways']:
				for ec in holder['ecs']:
					if not pathwayname in cache:
						cache[pathwayname] = {}
					if not ec in cache[pathwayname]:
						cache[pathwayname][ec] = {}
						
					#if this is printing the results of 'solve', then get flux activities...
					#reactionID	name	rev	pathway	ec	equation	confidence	gpr
					if showfluxvalues and self.REACTION2FLUXVALUE.get(reaction, '.') != '0':
						cache[pathwayname][ec][ ('\t').join((reaction, name, str(reversible), pathwayname, ec, self.REACTION2FLUXVALUE.get(reaction, '.'), reactionequation)) ] = 1
					if not showfluxvalues:
						cache[pathwayname][ec][ ('\t').join((reaction, name, str(reversible), pathwayname, ec, reactionequation, prr, grr, ref)) ] = 1
		
		paths = cache.keys()
		paths.sort()		
		
		if out:
			outfi = open(out, 'w')
				
			for path in paths:
				ecs = cache[path].keys()
				ecs.sort()
				for ec in ecs:
					for r in cache[path][ec]:
						print >>outfi, r
				#print >>outfi, '\n'
				
		else:
			for path in paths:
				ecs = cache[path].keys()
				ecs.sort()
				for ec in ecs:
					for r in cache[path][ec]:
						print r
				#print '\n'
				
		
	def write_constraints (self, outfilename):
		"Write a pickled object containing current model reaction constraints. Specify the filename."
		pickle.dump(self.CONSTRAINTS, open(outfilename, 'wb'), -1)
		
		
	def load_constraints (self, infilename):
		"Load model constraints from a pickled object. Specify the filename. For a given reaction, overwrites any existing constraints if there is a constraint in the pickled object."
		constraints_holder = pickle.load(open(infilename, 'rb'))
		for constraint in constraints_holder:
			(lbound, ubound) = constraints_holder[constraint]
			self.CONSTRAINTS[constraint] = (lbound, ubound)
			
										
	def set_sources (self, sourcelist):
		"Add source fluxes for all sources."
		
		#assert len(sourcelist) > 0, 'Must specify 1 or more sources'
		self.SOURCES = sourcelist
		
		for species in sourcelist:
			rev = False
			notes = {'SUBSYSTEM: SourceFlux':1, 'EC: .':1}
			boundaryflux_ID = 'R_SRC_' + species[2:]
			boundaryspecies_ID = species[:-1] + 'b'
			eq = [[(boundaryspecies_ID, '1')], [(species, '1')]]
			
			name, compartment, charge, boundaryCondition = '.', '.', '.', 'true'
			
			cb.add_reaction(self, boundaryflux_ID, name + ' source flux', rev, notes, eq)
			#if species[-1] == 'e':
			#	cb.set_constraint(self, boundaryflux_ID, '0', str(1000.0 * float(self.VMAX)) )


	def set_escapes (self, escapelist):
		"Add escape fluxes for all escapes (if escapes=[], add escapes for all species in model)."
		
		self.ESCAPES = escapelist
		
		for species in escapelist:
			#do not add escapes for boundary metabolites (avoid infinite regression of escapes)
			if not species[-1] == 'b':
				rev = False
				notes = {'SUBSYSTEM: EscapeFlux':1, 'EC: .':1}
				boundaryflux_ID = 'R_ESC_' + species[2:]
				boundaryspecies_ID = species[:-1] + 'b'
				eq = [[(species, '1')], [(boundaryspecies_ID, '1')]]
				specdict = self.SPECIES[species]

				name, compartment, charge, boundaryCondition = '.', '.', '.', 'true'
				
				cb.add_reaction(self, boundaryflux_ID, name + ' escape flux', rev, notes, eq)
				#if species[-1] == 'e':
				#	cb.set_constraint(self, boundaryflux_ID, '0', str(1000.0 * float(self.VMAX)) )
	
	
	def set_exchanges (self, exchangelist):
		"Add exchange fluxes for all exchanges."
		
		self.EXCHANGES = exchangelist
		
		for species, lb, ub in exchangelist:
			#do not add escapes for boundary metabolites (avoid infinite regression of escapes)
			if not species[-1] == 'b':
				rev = True
				notes = {'SUBSYSTEM: ExchangeFlux':1, 'EC: .':1}
				boundaryflux_ID = 'R_EXCH_' + species[2:]
				boundaryspecies_ID = species[:-1] + 'b'
				eq = [[(species, '1')], [(boundaryspecies_ID, '1')]]
				specdict = self.SPECIES[species]

				name, compartment, charge, boundaryCondition = '.', '.', '.', 'true'
				
				cb.add_reaction(self, boundaryflux_ID, name + ' exchange flux', rev, notes, eq)
				cb.set_constraint(self, boundaryflux_ID, lb, ub)
	
									
	def build(self, model_file, readquiet):
		#read and build initial model (just the reactions specified)
		if not readquiet:
			print 'model from', model_file
		file = open(model_file)
		while True:
			line = file.readline()
			if line == '': break
			line = line.rstrip()
			if line == '': continue
			if line[0] == '#': continue
			if line[:2] == 'R_' or line[0] == 'R':
				col = line.split('\t')
				if 'R_ILL_' in col[0]:
					continue
				
				[id, name, rev, pathwaysstr, ecsstr, stringequation] = col[0:6]
				
				reversibility, equation = eq_current.parse(stringequation)
				
				notes, pathways, ecs = {}, pathwaysstr.split('; '), ecsstr.split('; ')
				#if len(pathways) > 1: print id, "associated with > 1 pathways; splitting list on '; '"
				for pathway in pathways:
					notes['SUBSYSTEM: ' + pathway] = 1
				#if len(ecs) > 1: print id, "associated with > 1 ec numbers; splitting list on '; '"
				for ec in ecs:
					notes['EC: ' + ec] = 1
				
				#check read in of model...		
				#print stringequation
				#print '  ', id, name, reversibility, notes, equation
				
				cb.add_reaction(self, id, name, reversibility, notes, equation)
			#insert new code here to handle gpr, notes, refs, etc...
		if not readquiet:
			print


	def biomass(self, biomass_file, readquiet):
		#read and define biomass equation...
		if not readquiet:
			print 'biomass from', biomass_file
		biomass_equation = [[], []]
		file = open(biomass_file)
		while True:
			line = file.readline()
			if line == '': break
			line = line.rstrip()
			if line == '': continue
			if line[0] == '#': continue
			col = line.split('\t')
			[rawmet, coef, side] = col[0:3]
			
			#convert metabolites from forms like 'leu-L[c]' to 'M_leu_DASH_L_c'...
			met = eq_current.convert_metabolite_ext2int(rawmet)
			
			if 'reactant' in side:
				biomass_equation[0].append((met, coef))
				self.NOTSOURCES.append(met)
			elif 'product' in side:
				biomass_equation[1].append((met, coef))
				self.NOTSOURCES.append(met)
			else:
				print 'Warning: biomass component skipped because cannot be idenified as reactant or product-->', line
		#add biomass equation to model, and set this as the objective...
		cb.add_reaction(self, 'R_biomass_target', 'BiomassRxn', False, {'SUBSYSTEM: BiomassObjective':1, 'EC: .':1}, biomass_equation)
		cb.set_objective(self, 'Maximize', 'R_biomass_target')
		#cb.set_constraint(self, 'R_biomass_target', '0.0', str(float(self.VMAX)) )
		if not readquiet:
			print eq_current.makestring(biomass_equation, False)
			print


	def sources(self, sources_file, readquiet):
		#read and define source metabolites.	
		if not readquiet:
			print 'sources from', sources_file
		sources = {}
		file = open(sources_file)
		while True:
			line = file.readline()
			if line == '': break
			line = line.rstrip()
			if line == '': continue
			if line[0] == '#': continue
			col = line.split('\t')
			rawmet = col[0]
			
			if not readquiet:
				print ' ', rawmet

			#convert metabolites from forms like 'leu-L[c]' to 'M_leu_DASH_L_c'...
			met = eq_current.convert_metabolite_ext2int(rawmet)

			sources[met] = 1													
		sourcelist = sources.keys()
		cb.set_sources(self, sourcelist)
		if not readquiet:
			print


	def escapes(self, escapes_file, readquiet):
		#read and define escape metabolites.
		if not readquiet:
			print 'escapes from', escapes_file	
		escapes = {}
		file = open(escapes_file)
		while True:
			line = file.readline()
			if line == '': break
			line = line.rstrip()
			if line == '': continue
			if line[0] == '#': continue
			col = line.split('\t')
			rawmet = col[0]
			
			if not readquiet:
				print ' ', rawmet

			#convert metabolites from forms like 'leu-L[c]' to 'M_leu_DASH_L_c'...
			met = eq_current.convert_metabolite_ext2int(rawmet)

			escapes[met] = 1													
		escapelist = escapes.keys()
		cb.set_escapes(self, escapelist)
		if not readquiet:
			print


	def exchanges(self, exchanges_file, readquiet):
		#read and define exchange metabolites.
		if not readquiet:
			print 'exchanges from', exchanges_file	
		exchanges = {}
		file = open(exchanges_file)
		while True:
			line = file.readline()
			if line == '': break
			line = line.rstrip()
			if line == '': continue
			if line[0] == '#': continue
			col = line.split()					## starting here, some changes here to cope with MM2 output
			try:
				lb, ub = col[1], col[2]
			except:
				lb, ub = '-1000', '1000'
			rawmet = col[0]
			
			if not readquiet:
				print ' ', rawmet

			#convert metabolites from forms like 'leu-L[c]' to 'M_leu_DASH_L_c'...
			met = eq_current.convert_metabolite_ext2int(rawmet)

			exchanges[(met, lb, ub)] = 1													
		exchangelist = exchanges.keys()
		cb.set_exchanges(self, exchangelist)
		if not readquiet:
			print


	def constraints(self, constraints_file, readquiet):
		#read and define user-specified reaction constraints
		if not readquiet:
			print 'constraints from', constraints_file
		file = open(constraints_file)
		while True:
			line = file.readline()
			if line == '': break
			line = line.rstrip()
			if line == '': continue
			if line[0] == '#': continue
			col = line.split('\t')
			rxn, lbound, ubound = col[0], col[1], col[2]
			if not readquiet:
				print ' ', lbound, '<->', ubound, '\t', rxn
			cb.set_constraint(self, rxn, lbound, ubound)
		if not readquiet:
			print


	def notes(self, notes_file, readquiet):
		#read and define user-specified reaction constraints
		if not readquiet:
			print 'notes from', notes_file
		file = open(notes_file)
		pmidsearch = ''
		while True:
			line = file.readline()
			if line == '': break
			line = line.rstrip()
			if line == '': continue
			if line[0] == '#': continue
			col = line.split('\t')
			rxn, note = col[0], col[1]
			if 'PMID: ' == note[:6]:
				pmid = note.split()[1][:-1]
				pmidsearch = pmidsearch + pmid + ' '
			cb.add_note(self, rxn, note)
		#if not pmidsearch == '':
		#	print 'PMIDs:', pmidsearch
		if not readquiet:
			print 


	def gpr2(self, filename, readquiet):
		#read gpr.txt file
		if not readquiet:
			print 'gpr from', filename
		file = open(filename)
		while True:
			line = file.readline()
			if line == '': break
			line = line.rstrip()
			if line == '': continue
			if 'rg\t' == line[:3]:
				rg, rxn, gpr = line.split('\t')[0], line.split('\t')[1], line.split('\t')[2]
				cb.add_note(self, rxn, 'Gene_association: ' + gpr)
				#skip if there is no gpr...
				if gpr == '.':
					continue
				#if there is a gpr statement, add to REACTS and SIMPLEGPR...
				self.REACTS[rxn] = 1
				self.SIMPLEGPR[rxn] = gpr
				#if there is an 'and' statement in gpr, then this must be a protein complex...
				if 'and' in gpr:
					self.COMPLEXES[rxn] = 1
				#if there is more than one gene in gpr, and only 'or' statements, then this must be an isozyme...
				if len(gpr.split()) > 1 and not 'and' in gpr:
					self.ISOZYMES[rxn] = 1
				for item in gpr.split(' '):
					if item == 'and' or item == 'or':
						continue
					if item[0] == '(':
						item = item[1:]
					if item[-1] == ')':
						item = item[:-1]
					self.GENES[item] = 1
		if not readquiet:
			print


	def build_from_textfiles(self, modelfile, biomassfile=None, sourcesfile=None, escapesfile=None, exchangesfile=None, constraintsfile=None, notesfile=None, gprfile=None, readquiet=False):
		#one line command to build model from text files.
		cb.build(self, modelfile, readquiet)
		if biomassfile:
			cb.biomass(self, biomassfile, readquiet)
		if sourcesfile:
			cb.sources(self, sourcesfile, readquiet)
		if escapesfile:
			cb.escapes(self, escapesfile, readquiet)
		if exchangesfile:
			cb.exchanges(self, exchangesfile, readquiet)
		if constraintsfile:
			cb.constraints(self, constraintsfile, readquiet)
		if notesfile:
			cb.notes(self, notesfile, readquiet)
		if gprfile:
			cb.gpr2(self, gprfile, readquiet)



	def build_from_mm2(self, mm2file, readquiet=False):
		#this is command to build model from modelfile downloaded from mm2 (includes exchanges, gpr, model)
		modelfilename = mm2file[:-9] + '.model.txt'
		modelfile_ = open(modelfilename, 'w')
		gprfilename = mm2file[:-9] + '.gpr.txt'
		gprfile_ = open(gprfilename, 'w')
		exchangesfilename = mm2file[:-9] + '.exchanges.txt'
		exchangesfile_ = open(exchangesfilename, 'w')
		file = open(mm2file)
		while True:
			line = file.readline()
			if line == '': break
			line = line.rstrip()
			if line == '': continue
			if line[0] == '#': continue
			if line[:1] == 'R':
				print >>modelfile_, line
				continue
			elif line[:2] == 'rg':
				print >>gprfile_, line
				continue
			else:
				print >>exchangesfile_, line
				continue
		modelfile_.close()
		exchangesfile_.close()
		gprfile_.close()
		cb.build_from_textfiles(self, modelfilename, exchangesfile=exchangesfilename, gprfile=gprfilename, readquiet=readquiet)
		


		
	#calculator: given vector of gene presence/absence, calculate reaction presence/absence
	def calc (self):
		#make copies of the transcr, prots, reacts dictionaries so you don't change the global ones
		genes = self.GENES.copy()
		transcr = self.TRANSCR.copy()
		for t in transcr:
			#here and below, ensure that transcr[t] (or prots[p], ...) is a string so can use eval()
			#might be an integer, if you've specified a vector for presence/absence of transcripts
			transcr[t] = str(transcr[t])
			transcr[t] = eval(transcr[t])
		prots = self.PROTS.copy()
		for p in prots:
			prots[p] = str(prots[p])
			prots[p] = eval(prots[p])
		reacts = self.REACTS.copy()
		deletedrxns = {}
		for r in reacts:
			reacts[r] = str(reacts[r])
			reacts[r] = eval(reacts[r])
			if reacts[r] == 0 and r in self.REACTIONS:
				deletedrxns[r] = 1
		return deletedrxns
				

	def deletions (self, level):
		"Delete all genes, proteins, or reactions, one at a time."
		cb.solve(self, verbose=False)
		fullstatus, fullobjectivevalue = self.STATUS, self.OBJECTIVE_VALUE
		#item is a gene, protein, or reaction; level is self.GENES, self.PROTS, or self.REACTS
		lethals = {}
		for item in level:
			#delete item
			level[item] = 0
			#calculate consequences of deleted item according to boolean rules (i.e., which reactions are eliminated?)
			deletedrxns = cb.calc(self)
			#now constrain each reaction that is deleted by the change to have zero flux, attempt fba
			for r in deletedrxns:
				name, rev, notes, equation = self.REACTIONS[r]
				if rev:
					lbound = '-' + self.VMAX
				else:
					lbound = '0'
				default_lbound, default_ubound = self.CONSTRAINTS.get(r, (lbound, self.VMAX))
				cb.set_constraint(self, r, 0, 0)
				cb.solve(self, verbose=False)
				#if status != OPTIMAL or objective value is < 25% of 'wild type', print item, reaction, and results
				if (not self.STATUS == 'OPTIMAL') or (float(self.OBJECTIVE_VALUE) < 0.25 * float(fullobjectivevalue)):
					reactionequation = eq_current.makestring(equation, rev)
					print item + '\t' + r + '\t' + reactionequation + '\t' + self.STATUS + '\t' + str(self.OBJECTIVE_VALUE)
					lethals[item] = 1
				#reset reaction constraints to default	
				cb.set_constraint(self, r, default_lbound, default_ubound)
			#make item (gene, protein, ...) available again
			level[item] = 1
		return lethals
			

	def paul (self):
		"Algorithm to find a minimal set of source / escape reactions to convert a model with no solutions to one with a solution."
		#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		#WRITE *.dat FILE....

		datfilename = self.MODEL_ID + '.dat'
		datfile = open(datfilename, 'w')
		print >>datfile, "set REACTIONS := "
		for r in self.REACTIONS:
			print >>datfile, r,
		print >>datfile, ";\n"
		print >>datfile, "set METABOLITES := "
		for s in self.SPECIES:
			print >>datfile, s,
		print >>datfile, ";\n"
		print >>datfile, "set KNOWNSOURCES := "
		for s in self.SOURCES:
			print >>datfile, s[:-1] + 'b',
		print >>datfile, ";\n"
		print >>datfile, "set NOTSOURCES := "
		for n in self.NOTSOURCES:
			print >>datfile, n,
		print >>datfile, ";\n"
		print >>datfile, "set KNOWNESCAPES := "
		for e in self.ESCAPES:
			print >>datfile, e[:-1] + 'b',
		print >>datfile, ";\n"
		print >>datfile, "set NOTESCAPES := "
		for n in self.NOTESCAPES:
			print >>datfile, n,
		print >>datfile, ";\n"
		print >>datfile, "set BIOMASSREACTION := " + self.OBJECTIVE[1] + ";"
		print >>datfile, "\n"
		print >>datfile, "param S := "
		
		mets = {}
		constraints = {}
		for ID in self.REACTIONS:
					
			name, reversible, notes, equation = self.REACTIONS[ID]
						
			#if there have been specific constraints already set, use them...
			if ID in self.CONSTRAINTS:
				lbound, ubound = self.CONSTRAINTS[ID]
			#otherwise create constraints on the fly, using reaction reversibility and self.VMAX...
			else:
				ubound = self.VMAX
				if bool(reversible):
					lbound = '-' + self.VMAX
				else:
					lbound = '0'
			constraints[ID] = (lbound, ubound)	
				
			#create a data structure called 'mets': keys are metabolites, values are (reactionID, coef); this is essentially "S * v"
			for reactant in equation[0]:
				species, coef = reactant[0], '-' + str(reactant[1])
				#if '_b' == species[-2:]:
				#	continue
				if not species in mets:
					mets[species] = []
				mets[species].append((ID, coef))
				
			for product in equation[1]:
				species, coef = product[0], product[1]
				#if '_b' == species[-2:]:
				#	continue
				if not species in mets:
					mets[species] = []
				mets[species].append((ID, coef))

		#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

		for m in mets:
			rclist = mets[m]
			for rxn, coef in rclist:
				print >>datfile, m, rxn, coef
		print >>datfile, ";"
		print >>datfile, "param lb := "
		for r in self.REACTIONS:
			if r in constraints:
				lbound, ubound = constraints[r]
				print >>datfile, r, lbound
		print >>datfile, ";"
		print >>datfile, "param ub := "
		for r in self.REACTIONS:
			if r in constraints:
				lbound, ubound = constraints[r]
				print >>datfile, r, ubound
		print >>datfile, ";"
		print >>datfile, "param minbiomass := " + '1' + ";"
		print >>datfile, "end;\n"
		
		datfile.close()
		
		com1 = "/Users/seth/.lib/python/glpsol --tmlim 50 -m /Users/seth/.lib/python/metmodel/transport.mod -d " + datfilename + " -o transport.out"
		com2 = "rm " + datfilename
		com3 = "python /Users/seth/.lib/python/metmodel/out2transports.py transport.out" # > " + self.MODEL_ID + ".se"
		com4 = "rm transport.out"

		os.system(com1)
		#os.system(com2)
		os.system(com3)
		#os.system(com4)
		
	
		
	def ddeletions (self):
		"Double deletions at the reaction level."
		cb.solve(self, verbose=False)
		essential_reactions = cb.deletion_testing(self)

		fullstatus, fullobjectivevalue = self.STATUS, self.OBJECTIVE_VALUE
		for i, r in enumerate(self.REACTIONS.keys()):
			if not r == self.OBJECTIVE[1] and not r in essential_reactions and not 'R_ESC' in r and not 'R_SRC' in r:
				name, rev, notes, equation = self.REACTIONS[r]
				if rev:
					lbound = '-' + self.VMAX
				else:
					lbound = '0'
				default_lbound, default_ubound = self.CONSTRAINTS.get(r, (lbound, self.VMAX))
				cb.set_constraint(self, r, 0, 0)
				
				for j, r2 in enumerate(self.REACTIONS.keys()[i+1:]):
					if not r2 == self.OBJECTIVE[1] and not r2 in essential_reactions and not 'R_ESC' in r2 and not 'R_SRC' in r2:
						name2, rev2, notes2, equation2 = self.REACTIONS[r2]
						if rev2:
							lbound2 = '-' + self.VMAX
						else:
							lbound2 = '0'
						default_lbound2, default_ubound2 = self.CONSTRAINTS.get(r2, (lbound2, self.VMAX))
						cb.set_constraint(self, r2, 0, 0)
				
						cb.solve(self, verbose=False)
						
						if (not self.STATUS == 'OPTIMAL') or (float(self.OBJECTIVE_VALUE) < 0.25 * float(fullobjectivevalue)):
							reactionequation1 = eq_current.makestring(equation, rev)
							reactionequation2 = eq_current.makestring(equation2, rev2)
							print self.STATUS + '\t' + str(self.OBJECTIVE_VALUE)
							print r + '\t' + reactionequation1
							print r2 + '\t' + reactionequation2
							print
						cb.set_constraint(self, r2, default_lbound2, default_ubound2)
						
				cb.set_constraint(self, r, default_lbound, default_ubound)	

			
			
		
			
