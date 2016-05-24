is.generic.metacyc <-
function(equation, compound_df, id_col='UNIQUE.ID', smiles_col='SMILES', parent_col='TYPES', direction_type = c(' <=> ', ' => ')) {
  
  # Check parameters
  if(id_col %in% colnames(compound_df) == FALSE) {
    stop(sprintf('There is no column called "%s" for id_col. Please provide proper parameter.', id_col))
  }
  if(smiles_col %in% colnames(compound_df) == FALSE) {
    stop(sprintf('There is no column called "%s" for smiles_col. Please provide proper parameter.', smiles_col))
  } 
  
  result = logical(length(equation))
  # We don't need localization for instantiation
  pattern_localization = '\\(.+\\)'
  testEquation = gsub(pattern_localization, '', equation)
  
  ind_n = grep('n ', testEquation)
  if(length(ind_n) >0) {
    message(sprintf('%s coefficient "n" containing equation(s) are found being processed by substituting 1-10 for n', length(ind_n)))
    tmpEquation = gsub('\\(n\\+2\\)', 3, equation[ind_n])
    tmpEquation = gsub('7n', 7, tmpEquation)
    tmpEquation = gsub('2n', 2, tmpEquation)
    tmpEquation = gsub('n ', '', tmpEquation)
    equation[ind_n] = tmpEquation
  }
  
  message('Build parent table')
  parentTable = build.subtable(compound_df, id_col, parent_col, '///')
  parentTable[is.na(parentTable)] = ''
  colnames(parentTable) = c("child", 'parent')
  
  # Remove too broad child-parent set
  ## 'Compounds'
  parentTable2 = parentTable[parentTable$parent != 'Compounds',]
  
  for(i in 1:length(equation)) {
    # Get Participant
    participants = .get.participant(equation[i], direction_type)
    participants = gsub('\\|','', participants)
    
    # Get compound class (in this version, only considered alkyl and partly polymerized compound)
    alkyl = compound_df[compound_df[[id_col]] %in% participants,][grep('R', compound_df[compound_df[[id_col]] %in% participants, smiles_col]),id_col]
    
    # Get compound class (in this version, only considered alkyl and partly polymerized compound)
    classCompound = participants[participants %in% parentTable2$parent]
    alkyl_or_class = unique(c(alkyl, classCompound))
    
    if(length(alkyl_or_class) > 0) {
      result[i] = TRUE
    }
  }
  return(result)
}
