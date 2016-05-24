ApplyGrammarRule <- function(current.codon, startSymb, grammar) {
  # applies current codon to the startSymb (which can be anystring) 
  # returns new string
  
  TERMINAL = TRUE
  NON_TERMINAL = FALSE
  
  sep.symbs = GetFirstNonTerminalandRest(startSymb)
  if (sep.symbs[[1]] == NON_TERMINAL){
    return (startSymb) # no more replacements
  }
  
  terminal.symb = sep.symbs[[1]]
  non.terminal.symb = sep.symbs[[2]]
  rest = sep.symbs[[3]]

  possible.choices = GetPossibleRuleChoices(non.terminal.symb, grammar)
  choice.no = current.codon %% possible.choices
  chosen.rule = ChosenGrammarRule(non.terminal.symb, choice.no+1, grammar)

  new.startSymb = paste0(terminal.symb, chosen.rule, rest)
  new.startSymb
}

ChosenGrammarRule <- function(rule.name, choice, grammar.rules){
  rule.index = GetGrammarRuleIndex(grammar.rules, rule.name)
  rule.to.apply = grammar.rules$def[[rule.index]][[2]][[choice]]
  return(rule.to.apply)
}

GetPossibleRuleChoices <- function(non.terminal.symb, grammar.rules){
  if (!is.character(non.terminal.symb))
    stop("Invalid non-string input to GetPossibleRuleChoices()")
  ruleInd = GetGrammarRuleIndex(grammar.rules, paste(non.terminal.symb))
  return(grammar.rules$ruleSizes[ruleInd])
}

GetFirstNonTerminalandRest <- function(rule){
  
  NON_TERMINAL = FALSE
  
  grep.res.start = gregexpr("<", rule)[[1]]
  if (grep.res.start[1] == -1){
    return(list(NON_TERMINAL))
  }

  grep.res.end = gregexpr(">", rule)[[1]]
  if (grep.res.end[1] == -1){
    stop("invalid rule representation")
  }

  terminals = substr(rule, 1, grep.res.start-1)
  first.nonterminal = substr(rule, grep.res.start+1, grep.res.end-1)
  rest = substr(rule, grep.res.end+1, nchar(rule))

  return(list(terminals = terminals, 
              nonterminal = first.nonterminal, 
              rest = rest))
}

GetPossibleRuleChoicesFirstSymbol <- function (symbol, grammar) {
  
  TERMINAL = TRUE
  NON_TERMINAL = FALSE
  
  sep.symbs = GetFirstNonTerminalandRest(symbol)
  if (sep.symbs[[1]] == NON_TERMINAL){
    return (0) # no more replacements
  }
  
  possible.choices = GetPossibleRuleChoices(sep.symbs$nonterminal, grammar)
  
  return (possible.choices)
}

IsSymbolTerminal <- function(symb){
  # checks if the grammar is terminal or not
  grep.res.start = gregexpr("<", symb)[[1]]
  return (grep.res.start[1] == -1)
}

IsGrammarTerminal <- IsSymbolTerminal

GetGrammarRuleIndex <- function(grammar, rule.name) {
  # finds the index of rule in the grammar
  attr(grammar$defIndex, rule.name)
}

TraverseCodon <- function(codon.string, startSymb, grammar) {
  
  total.genes = length(codon.string)
  
  for (i in 1:total.genes){
    res = ApplyGrammarRule(codon.string[i], startSymb, grammar)
    if (res == startSymb){ # if not replacement is performed, we have reached a terminal
      return(startSymb) 
    }
    startSymb = res  
  }
  startSymb
}

GetNumNonTerminals <- function(startSymb) {
  x = 0
  while (!IsSymbolTerminal(startSymb)) {
    x = x + 1
    sep.symb = GetFirstNonTerminalandRest(startSymb)    
    startSymb = sep.symb$rest
  }
  
  return(x)
}
