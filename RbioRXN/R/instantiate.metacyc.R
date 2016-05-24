instantiate.metacyc <-
function(equation, chemical_table, id_col='UNIQUE.ID', parent_col='TYPES', formula_col='CHEMICAL.FORMULA', smiles_col='SMILES', inchi_col='INCHI', direction_type=c(" => ", " <=> ")) {
	result = list()
  
	# Substitute 1 for n
	pattern_localization = '\\(.+\\)'
	testEquation = gsub(pattern_localization, '', equation)

	ind_n = grep('n ', testEquation)
	if(length(ind_n) >0) {
		cat(sprintf('%s coefficient "n" containing equations are found being processed by substituting 1-10 for n\n', length(ind_n)))
		for(n in 1:10) {
			tmpEquation = gsub('\\(n\\+2\\) ', paste(n+2, ' ', sep=''), equation[ind_n])
			tmpEquation = gsub('7n ', paste(7*n, ' ', sep=''), tmpEquation)
			tmpEquation = gsub('2n ', paste(2*n, ' ', sep=''), tmpEquation)
			tmpEquation = gsub('n ', paste(n, ' ', sep=''), tmpEquation)
			equation = c(equation, tmpEquation)
		}
		equation = equation[-ind_n]
	}
	# Get Participant
	participants = .get.participant(equation, direction_type=direction_type)
  participants2 = gsub('\\|','', participants)
  
  # Build parent-child table
	cat('Build parent table\n')
	parentTable = build.subtable(chemical_table, id_col, parent_col, '///')
	parentTable[is.na(parentTable)] = ''
	colnames(parentTable) = c("child", 'parent')
  
	# Remove too broad child-parent set
  ## 'Compounds'
	parentTable2 = parentTable[parentTable$parent != 'Compounds',]
  
	# Get compound class (in this version, only considered alkyl and partly polymerized compound)
	classCompound = participants2[participants2 %in% parentTable$parent]
  
	if(length(classCompound) == 0) {
		stop('There is no class compound in your equations')
	}

	cat(sprintf('%i class compounds found\n', length(classCompound)))

	# Instantiation
	## Determine direction symbol
	directionalities = numeric(length(equation))
	for(i in direction_type) {
		ind_direction = grepl(i, equation) 
		directionalities[ind_direction] = i
	}
  
	number_of_instantiated_generic = 0
  for(i in 1:length(equation)) {
    instantiated_reactions = tryCatch({
      instantiate(equation[i], chemical_table, parentTable2, directionalities[i], classCompound, id_col, inchi_col, smiles_col, formula_col)
    }, error = function(cond) {
      cat('ERROR:\n')
      print(cond)
      return(0)
    })
    if(mode(instantiated_reactions) == 'character') {
      number_of_instantiated_generic = number_of_instantiated_generic + 1
    }
    result[[equation[i]]] = instantiated_reactions
  }
	cat(sprintf('# of instantiated generic reaction: %i\n', number_of_instantiated_generic))
	return(result)
}
