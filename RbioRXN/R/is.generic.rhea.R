is.generic.rhea <-
function(equation, chebi_df, id_col='chebi', parent_col='parent', formula_col='formula', direction_type = c(' <=> ', ' => ', ' <\\?> ')) {
  
  # Check parameters
  if(id_col %in% colnames(chebi_df) == FALSE) {
    stop(sprintf('There is no column called "%s" for id_col. Please provide proper parameter.', id_col))
  }
  if(parent_col %in% colnames(chebi_df) == FALSE) {
    stop(sprintf('There is no column called "%s" for parent_col. Please provide proper parameter.', parent_col))
  }
  if(formula_col %in% colnames(chebi_df) == FALSE) {
    stop(sprintf('There is no column called "%s" for formula_col. Please provide proper parameter.', formula_col))
  } 
  
  result = logical(length(equation))
  # We don't need localization for instantiation
  pattern_localization = '\\(.+\\)'
  testEquation = gsub(pattern_localization, '', equation)
  
  ind_n = grep('n', testEquation)
  if(length(ind_n) >0) {
    message(sprintf('%s coefficient "n" containing equations are found being processed by substituting 1-10 for n', length(ind_n)))
    tmpEquation = gsub('\\(n\\+2\\)', 3, equation[ind_n])
    tmpEquation = gsub('7n', 7, tmpEquation)
    tmpEquation = gsub('2n', 2, tmpEquation)
    tmpEquation = gsub('n ', '', tmpEquation)
    equation[ind_n] = tmpEquation
  }
  
  for(i in 1:length(equation)) {
    # Get Participant
    participants = .get.participant(equation[i], direction_type)
    
    # Get compound class (in this version, only considered alkyl and partly polymerized compound)
    alkyl = chebi_df[chebi_df[[id_col]] %in% participants,][grep('R', chebi_df[chebi_df[[id_col]] %in% participants, formula_col]),id_col]
    partly_polymerized = chebi_df[chebi_df[[id_col]] %in% participants,][grep('\\)n', chebi_df[chebi_df[[id_col]] %in% participants, formula_col]),id_col]
    classCompound = unique(c(alkyl, partly_polymerized))
    
    if(length(classCompound) > 0) {
      result[i] = TRUE
    }
  }
  return(result)
}
