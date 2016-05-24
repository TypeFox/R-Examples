GrammarNumOfExpressions <- 
  function(grammar, 
           max.depth = GrammarGetDepth(grammar), 
           startSymb = GrammarStartSymbol(grammar),
           ...) {
  # returns the maximum number of expressions generatable by grammar
  # max.depth is the maximum depth to search into the tree. It is used to avoid
  # looping in recursive grammar
  
  if ((!IsSymbolTerminal(startSymb)) & (max.depth == 0))
    return (0)
  
  depth = 1
  while(IsSymbolTerminal(startSymb) == FALSE) {
    sep.symbs = GetFirstNonTerminalandRest(startSymb)
    possible.choices = GetPossibleRuleChoices(sep.symbs$nonterminal, grammar)
    
    depth_list = sapply(1:possible.choices, function(choice.no) {
      chosen.rule = ChosenGrammarRule(sep.symbs$nonterminal, choice.no, grammar)
      GrammarNumOfExpressions(grammar, max.depth - 1, startSymb = chosen.rule)
    })
    depth = depth * sum(depth_list)
    
    startSymb = sep.symbs$rest
  }
  
  return (depth)
}

GetGrammarNumOfExpressions <- GrammarNumOfExpressions
