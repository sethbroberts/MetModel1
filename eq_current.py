#script purpose: reaction equation parser and printer
"""
Parse reactions of the following form:
[c] : akg + asp-L <==> glu-L + oaa
[c] : 2cpr5p + h --> 3ig3p + co2 + h2o
[c] : 2 ala-D + atp <==> adp + alaala + h + pi
atp[c] + h2o[c] + urea[e] --> adp[c] + h[c] + urea[c] + pi[c]
[c] : gtp + 3 h2o --> 25dhpp + for + 2 h + ppi

Into data structures like this:
True, [[('M_akg_c', '1'), ('M_asp_DASH_L_c', '1')], [('M_glu_DASH_L_c', '1'), ('M_oaa_c', '1')]]
False, [[('M_2cpr5p_c', '1'), ('M_h_c', '1')], [('M_3ig3p_c', '1'), ('M_co2_c', '1'), ('M_h2o_c', '1')]]
True, [[('M_ala_DASH_D_c', '2'), ('M_atp_c', '1')], [('M_adp_c', '1'), ('M_alaala_c', '1'), ('M_h_c', '1'), ('M_pi_c', '1')]]
False, [[('M_atp_c', '1'), ('M_h2o_c', '1'), ('M_urea_e', '1')], [('M_adp_c', '1'), ('M_h_c', '1'), ('M_urea_c', '1'), ('M_pi_c', '1')]]
False, [[('M_gtp_c', '1'), ('M_h2o_c', '3')], [('M_25dhpp_c', '1'), ('M_for_c', '1'), ('M_h_c', '2'), ('M_ppi_c', '1')]]

First member of tuple tells whether reaction is reversible, second looks like this:
[[(reactant1, coefficient_reactant1), (reactant2, coefficient_reactant2), ...], [(product1, coefficient_product1), (product2, coefficient_product2), ...]]
"""


def no_compartment(rxnequation, splitter):
	print rxnequation, splitter
	raw_equation_array = [rxnequation.split(splitter)[0].split(' + '), rxnequation.split(splitter)[1].split(' + ')]
	equation_array = [[], []]
	for i, side in enumerate(raw_equation_array):
		for speccoef in side:
			tmparray = speccoef.split()
			assert len(tmparray) < 3, "Missing a '+' sign? %s" % rxnequation
			if len(tmparray) == 1:		#only species name is in tmparray, implies stoichiometric coefficient is '1'
				species, coefficient = 'M_' + tmparray[0], '1'
			else:						#both stoichiometric coefficient and species name are in tmparray
				species, coefficient = 'M_' + tmparray[1], tmparray[0]
			
			equation_array[i].append((species, coefficient))
	if splitter == ' <-- ':		
		return [equation_array[1], equation_array[0]]
	else:
		return equation_array

def single_compartment(compartment_suffix, rxnequation, splitter):
	raw_equation_array = [rxnequation.split(splitter)[0].split(' + '), rxnequation.split(splitter)[1].split(' + ')]
	equation_array = [[], []]
	for i, side in enumerate(raw_equation_array):
		for speccoef in side:
			tmparray = speccoef.split()
			assert len(tmparray) < 3, "Missing a '+' sign? %s" % rxnequation
			if len(tmparray) == 1:		#only species name is in tmparray, implies stoichiometric coefficient is '1'
				species_raw, coefficient = 'M_' + tmparray[0] + compartment_suffix, '1'
			else:						#both stoichiometric coefficient and species name are in tmparray
				species_raw, coefficient = 'M_' + tmparray[1] + compartment_suffix, tmparray[0]
			species = species_raw.replace('-', '_DASH_')
			
			equation_array[i].append((species, coefficient))
	if splitter == ' <-- ':		
		return [equation_array[1], equation_array[0]]
	else:
		return equation_array


def multi_compartment(rxnequation, splitter):
	raw_equation_array = [rxnequation.split(splitter)[0].split(' + '), rxnequation.split(splitter)[1].split(' + ')]
	equation_array = [[], []]
	for i, side in enumerate(raw_equation_array):
		for speccoef in side:
			tmparray = speccoef.split()
			assert len(tmparray) < 3, "Missing a '+' sign? %s" % rxnequation
			if len(tmparray) == 1:		#only species name is in tmparray, implies stoichiometric coefficient is '1'
				species_raw, coefficient = 'M_' + tmparray[0][:-3] + '_' + tmparray[0][-2:-1], '1'
			else:						#both stoichiometric coefficient and species name are in tmparray
				species_raw, coefficient = 'M_' + tmparray[1][:-3] + '_' + tmparray[1][-2:-1], tmparray[0]
			species = species_raw.replace('-', '_DASH_')
			
			equation_array[i].append((species, coefficient))
	if splitter == ' <-- ':		
		return [equation_array[1], equation_array[0]]
	else:
		return equation_array


def parse (rxnequation):
	#two main decisions that will determine how rxnequation is parsed:
	#1. is it reversible or not?
	#2. is it within a single compartment or not?

	compartment_suffix, reaction_arrow = '', ''

	if ' <==> ' in rxnequation:
		#reversible: split rxnequation on ' <==> '
		reaction_arrow = ' <==> '
		reversibility = True
	elif ' --> ' in rxnequation:
		#irreversible: split rxnequation on ' --> '
		reaction_arrow = ' --> '
		reversibility = False
	elif ' <-- ' in rxnequation:
		#irreversible: split rxnequation on ' <-- '; following initial parsing, this will be written left to right
		reaction_arrow = ' <-- '
		reversibility = False
	elif ' <=> ' in rxnequation:
		#reversible: split rxnequation on ' <=> '; this arrow is in .wil files
		reaction_arrow = ' <=> '
		reversibility = True



	assert reaction_arrow in rxnequation, "Mistake in reaction direction arrow. %s" % rxnequation


	#are compartments included in equation?
	if '[' in rxnequation:
		if rxnequation[0] == '[':
			#reaction is entirely within one compartment
			rxnarray = rxnequation.split()
			compartment_suffix = '_' + rxnarray[0][1:2]
			rxnequation = (' ').join(rxnarray[2:])
			equation = single_compartment(compartment_suffix, rxnequation, reaction_arrow)	
		else:
			#reaction in >1 compartments
			equation = multi_compartment(rxnequation, reaction_arrow)
	else:
		equation = no_compartment(rxnequation, reaction_arrow)
	
	return reversibility, equation	
		
		
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


def makestring (eq, rev):
	#Pretty print the reaction equation
	
	#from reversibility, assign reaction arrow
	if bool(rev):
		arrow = '<==>'
	else:
		arrow = '-->'
	
	#determine whether reaction occurs in single or multiple compartments
	compartments = {}
	for side in eq:
		for spec, coef in side:
			compartments[spec[-1]] = 1
	if len(compartments.keys()) == 1:
		prefix = '[' + compartments.keys()[0] + '] : '
		oneCompartment = True
	else:
		prefix = ''
		oneCompartment = False
		
	#build string for reaction	
	stringed_eq = prefix
	for i, side in enumerate(eq):
		for spec, coef in side:
			if oneCompartment:
				spec = spec[2:-2]
			else:
				spec = spec[2:-2] + '[' + spec[-1] + ']'
			if '_DASH_' in spec:
				spec = spec.replace('_DASH_', '-')
			#take off useless trailing '.0' (i.e., '2.0' -> '2')
			if coef[-2:] == '.0':
				coef = coef[:-2]
			if coef == '1' or coef == '1.0':
				stringed_eq = stringed_eq + spec + ' + '
			else:
				stringed_eq = stringed_eq + coef + ' ' + spec + ' + '
		#if this is reactant side, place arrow after all reactants added to string
		if i == 0:
			stringed_eq = stringed_eq[:-3] + ' ' + arrow + ' '
	#chop off last ' + ' placed at end of equation
	stringed_eq = stringed_eq[:-3]

	return stringed_eq

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
def makeoldstring (eq, rev):
	#Pretty print the reaction equation in the 'old' style, e.g., 
	
	#from reversibility, assign reaction arrow
	if bool(rev):
		arrow = '<==>'
	else:
		arrow = '-->'
	
	#build string for reaction	
	stringed_eq = ''
	for i, side in enumerate(eq):
		for spec, coef in side:
			if coef == '1':
				stringed_eq = stringed_eq + spec + ' + '
			else:
				stringed_eq = stringed_eq + coef + ' ' + spec + ' + '
		#if this is reactant side, place arrow after all reactants added to string
		if i == 0:
			stringed_eq = stringed_eq[:-3] + ' ' + arrow + ' '
	#chop off last ' + ' placed at end of equation
	stringed_eq = stringed_eq[:-3]

	return stringed_eq

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def convert_metabolite_ext2int (rawmet):
	#convert metabolites from forms like 'leu-L[c]' to 'M_leu_DASH_L_c' (external to internal representation)...
	if '-' in rawmet:
		rawmet = rawmet.replace('-', '_DASH_')
	met = 'M_' + rawmet[:-3] + '_' + rawmet[-2]
	
	return met

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
def convert_metabolite_int2ext (rawmet):
	#convert metabolites from forms like 'M_leu_DASH_L_c' to 'leu-L[c]' (internal to external representation)...
	rawmet = rawmet.replace('_DASH_', '-')
	met = rawmet[2:-2] + '[' + rawmet[-1] + ']'
	
	return met

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
def makestring_nocomp (eq, rev):
	#Pretty print the reaction equation when there are no compartments 
	
	#from reversibility, assign reaction arrow
	if bool(rev):
		arrow = '<==>'
	else:
		arrow = '-->'
	
	#build string for reaction	
	stringed_eq = ''
	for i, side in enumerate(eq):
		for spec, coef in side:
			if coef == '1':
				stringed_eq = stringed_eq + spec[2:] + ' + '
			else:
				stringed_eq = stringed_eq + coef + ' ' + spec[2:] + ' + '
		#if this is reactant side, place arrow after all reactants added to string
		if i == 0:
			stringed_eq = stringed_eq[:-3] + ' ' + arrow + ' '
	#chop off last ' + ' placed at end of equation
	stringed_eq = stringed_eq[:-3]

	return stringed_eq


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def cac_eq_rewrite (rxnequation):

	transport = False
	#is this entirely in cytoplasm, or a transport reaction?
	if '(extracellular)' in rxnequation:
		transport = True
		
	reaction_arrow, reversibility = determine_reversibility(rxnequation)
	raw_equation_array = [rxnequation.split(reaction_arrow)[0].split(' + '), rxnequation.split(reaction_arrow)[1].split(' + ')]

	new_equation_arr = ['', '']

	for i, side in enumerate(raw_equation_array):
		for speccoef in side:
			tmparray = speccoef.split()
			#assert len(tmparray) < 3, "Missing a '+' sign? %s" % rxnequation
			if len(tmparray) == 1:		#only species name is in tmparray, implies stoichiometric coefficient is '1'
				coefficient = '1'
				species = tmparray[0]
			else:						#both stoichiometric coefficient and species name may be in tmparray
				#see if first element is a float; if it is, then this is coef, and rest is spec
				try:
					coefficient = str(float(tmparray[0]))
					species = (' ').join(tmparray[1:])
				#if first element is not a float, then coef is 1 and entire tmparray is spec
				except:
					coefficient = '1'
					species = speccoef
				#strip parens out of coef, and change 1.0 to 1...
				coefficient = coefficient.replace('(', '')
				coefficient = coefficient.replace(')', '')
				if coefficient == '1.0':
					coefficient = '1'
			
			species = species.lower()
			species = species.replace('[', '_LBRACK_')
			species = species.replace(']', '_LBRACK_')

			#modify species name based on whether or not this is a transport...
			if transport:
				if ' (extracellular)' in species:
					species = species.replace(' (extracellular)', '[e]')
				else:
					species = species + '[c]'

			#find and replace characters that glpk may choke on...
			species = species.replace(' ', '_')
			species = species.replace('(', '_LPARENS_')
			species = species.replace(')', '_RPARENS_')
			species = species.replace("'", '_APOST_')
			species = species.replace('-', '_DASH_')
			species = species.replace('+', '_PLUS_')

			if not coefficient == '1':
				new_equation_arr[i] = new_equation_arr[i] + coefficient + ' ' + species + ' + '
			else:
				new_equation_arr[i] = new_equation_arr[i] + species + ' + '
		
		new_equation_arr[i] = new_equation_arr[i][:-3]
	
	if reaction_arrow == ' <-- ':
		new_equation_arr = [new_equation_arr[1], new_equation_arr[0]]
		reaction_arrow = ' --> '

	new_equation = new_equation_arr[0] + reaction_arrow + new_equation_arr[1]
	
	if not transport:
		new_equation = '[c] : ' + new_equation
		
	return new_equation
	
	
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def determine_reversibility (rxnequation):
	reaction_arrow, reversibility = '.', '.'
	if ' <==> ' in rxnequation:
		#reversible: split rxnequation on ' <==> '
		reaction_arrow = ' <==> '
		reversibility = True
	elif ' --> ' in rxnequation:
		#irreversible: split rxnequation on ' --> '
		reaction_arrow = ' --> '
		reversibility = False
	elif ' <-- ' in rxnequation:
		#irreversible: split rxnequation on ' <-- '; following initial parsing, this will be written left to right
		reaction_arrow = ' <-- '
		reversibility = False
	elif ' = ' in rxnequation:
		#revisible: split on ' = '
		reaction_arrow = ' = '
		reversibility = True
	elif ' => ' in rxnequation:
		#irreversible: split on ' => '
		reaction_arrow = ' => '
		reversibility = False
	elif ' <=> ' in rxnequation:
		#reversible: split on ' <=> '
		reaction_arrow = ' <=> '
		reversibility = True
	elif ' <= ' in rxnequation:
		#irreversible: split on ' <= '
		reaction_arrow = ' <= '
		reversibility = False
	assert reaction_arrow in rxnequation, "Mistake in reaction direction arrow. %s" % rxnequation
	return reaction_arrow, reversibility
	
