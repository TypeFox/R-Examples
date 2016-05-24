CreateGrammar <- function(ruleDef, startSymb) {
  
  # if the ruledef is a file, read it
  if (("connection" %in% class(ruleDef)) | is.character(ruleDef)) {
    ruleDef = ReadBNFFile(ruleDef)
  } 
  # if the ruleDef is a symbolic expression
  if ("GERule" %in% class(ruleDef[[1]]))  {
    # convert it to list of text
    ruleDef = SymbolicRuleToListRule(ruleDef)
  }
  
  # trim brackets from the rule index
  for (i in seq_along(ruleDef)) {
    ruleDef[[i]][[1]] = trim_brackets(ruleDef[[i]][[1]])
  }
  
  # parse the text grammar
  grammar = CreateGrammarFromTextList(ruleDef)
  
  # add the starting symbol
  if (missing(startSymb)) {
    startSymb = as.character(grammar$defIndex[1])
  }
  grammar$startSymb = paste0('<', trim_brackets(startSymb), '>') # start Symbol
  
  class(grammar) = "grammar"
  
  return (grammar)
}

CreateGrammarFromTextList <- function(ruleDef) {
  # extract information about the grammar from list object
  # to create a grammar object
  
  # the name of expressions
  ruleDefIndex = sapply(ruleDef, function(r) r[[1]])
  
  for (i in 1:length(ruleDef)) # with their id
    attr(ruleDefIndex, ruleDef[[i]][[1]]) = i
  
  # the number of rules per expressions
  ruleSizes = sapply(ruleDef, function(r) length(r[[2]]))
  
  # attach
  grammar = list(def = ruleDef,             # Grammar in List Format
                 defIndex = ruleDefIndex,     # Index of each Rule
                 ruleSizes = ruleSizes,       # Length of each Rule as attr
                 maxRuleSize = max(ruleSizes))# maximum choice in a rule
  
  grammar
}

