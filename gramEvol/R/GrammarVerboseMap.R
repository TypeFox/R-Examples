GrammarVerboseMap <- function(inputString, grammar, wrappings = 3) {
  
  if (class(grammar) != "grammar") {
    stop("Invalid Grammar Class")
  }
  
  startSymb <- grammar$startSymb # start symbol

  # startubg step
  df_currentIteration = data.frame(
                  Step = 0,
                  Codon = NA,
                  Symbol = "",
                  Rule = "starting:",
                  Result = startSymb)
  
  df = NULL
  step = 1
  for (i in 1:wrappings) {

    if (i > 1) {
      cat("Wrapping string to position 0\n")
    }
    
    TERMINAL = TRUE
    NON_TERMINAL = FALSE
    
    codon.string = inputString
    
    total.genes = length(codon.string)
    
    for (j in 1:total.genes){
      
      sep.symbs = GetFirstNonTerminalandRest(startSymb)
      if (sep.symbs[[1]] == NON_TERMINAL){
        res = startSymb
      } else {
        terminal.symb = sep.symbs[[1]]
        non.terminal.symb = sep.symbs[[2]]
        rest = sep.symbs[[3]]
        
        possible.choices = GetPossibleRuleChoices(non.terminal.symb, grammar)
        choice.no = codon.string[j] %% possible.choices
        chosen.rule = ChosenGrammarRule(non.terminal.symb, choice.no+1, grammar)

        res = paste0(terminal.symb, chosen.rule, rest)
      }
      
      df_currentIteration = rbind(df_currentIteration, 
                     data.frame(Step = step,
                                Codon = codon.string[j],
                                Symbol = paste0('<', non.terminal.symb, '>'),
                                Rule = chosen.rule,
                                Result = res))
      if (res == startSymb){ # if not replacement is performed, we have reached a terminal
        break
      }
      startSymb = res
      step = step + 1
    }
    
    result.string = startSymb
    
    if (is.na(df_currentIteration[1,2])) {
      df_currentIteration[1,2] = ""
    }
    print(df_currentIteration, right=FALSE, row.names=FALSE)
    
    df = rbind(df, df_currentIteration)
    df_currentIteration = NULL

    if (IsSymbolTerminal(result.string)) {
      cat("Valid Expression Found\n")
      #cat(as.character(parse(text=unescape.gt.lt(result.string))), '\n')
      break
    } else {
      cat("Non-terminal expression\n")
      exprs = list(expr = result.string, type = "NT", parsed = NULL)
    }
  }  
  
  return (df)
}