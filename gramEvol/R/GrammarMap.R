GrammarMap <- function(inputString, 
                       grammar, 
                       wrappings = 3, 
                       verbose = FALSE) {
  ## Coverts an Integer String to Expression using Grammar
  
  if (class(grammar) != "grammar") {
    stop("Invalid Grammar Class")
  }
  if (!is.numeric(inputString)) {
    stop("Invalid Rule Sequence")
  }
  
  if (verbose) {
    GrammarVerboseMap(inputString, grammar, wrappings = 3)
  }
  
  result.string <- grammar$startSymb # start symbol
  
  for (i in 1:wrappings) {
    result.string = TraverseCodon(inputString, result.string, grammar)  
    
    if (IsSymbolTerminal(result.string)) {
      exprs = list(expr = result.string, type = "T", parsed = parse(text=unescape.gt.lt(result.string)))
      break
    } else {
      exprs = list(expr = result.string, type = "NT", parsed = NULL)
    }
  }  
  
  # also keep other information
  exprs$genotype = inputString
  exprs$wrappings = wrappings
  
  class(exprs) = "GEPhenotype"
  
  return (exprs)
}

GrammarGenotypeToPhenotype <- GrammarMap

