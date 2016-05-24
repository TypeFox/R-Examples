GrammarIsRecursive <- function(grammar, 
                               startSymb = GrammarStartSymbol(grammar),
                               ...) {

  # call the recursive function, start from no observed symbol
  RecursiveGrammarIsRecursive(grammar, startSymb, NULL)
}

RecursiveGrammarIsRecursive <- function(grammar, 
                               startSymb = GrammarStartSymbol(grammar), 
                               symb.list = NULL) {
  
  same.level.symbols = list()
  while (IsSymbolTerminal(startSymb) == FALSE) {
    sep.symbs = GetFirstNonTerminalandRest(startSymb)
    nonterminal = sep.symbs$nonterminal
    
    # check if the same symbol hasn't been processed in the same level
    if (!nonterminal %in% same.level.symbols) {
      
      # check if symbol is found before
      if (nonterminal %in% symb.list) {
        return (TRUE)
      }
      
      # otherwise, add to the seen list
      symb.list = c(symb.list, nonterminal)
      same.level.symbols = c(same.level.symbols, nonterminal)
      
      # dive into the symbol
      possible.choices = GetPossibleRuleChoices(nonterminal, grammar)
      is.recursive = sapply(1:possible.choices, function(choice.no) {
        chosen.rule = ChosenGrammarRule(nonterminal, choice.no, grammar)
        RecursiveGrammarIsRecursive(grammar, startSymb = chosen.rule, symb.list)
      })
      
      if (any(is.recursive)) {
        return (TRUE)
      }
    } 
    
    startSymb = sep.symbs$rest
  }
  
  return(FALSE)
}
