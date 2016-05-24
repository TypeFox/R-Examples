GrammarGetDepth <- function(grammar, 
                         max.depth = max(length(grammar$def), 4), 
                         startSymb = GrammarStartSymbol(grammar),
                         ...) {
  # returns the maximum depth of grammatical tree
  
  if (IsSymbolTerminal(startSymb) || (max.depth == 0))
    return (0)
  
  depth = 0
  while(IsSymbolTerminal(startSymb) == FALSE) {
    sep.symbs = GetFirstNonTerminalandRest(startSymb)
    possible.choices = GetPossibleRuleChoices(sep.symbs$nonterminal, grammar)
    
    depth_list = sapply(1:possible.choices, function(choice.no) {
      chosen.rule = ChosenGrammarRule(sep.symbs$nonterminal, choice.no, grammar)
      return (GrammarGetDepth(grammar, max.depth - 1, startSymb = chosen.rule))
    })
    depth = max(depth, depth_list)
    
    startSymb = sep.symbs$rest
  }
  
  return (1 + depth)
}

GetGrammarDepth <- GrammarGetDepth
