GrammarMaxSequenceRange <- function(grammar, 
                            max.depth = GrammarGetDepth(grammar), 
                            startSymb = GrammarStartSymbol(grammar),
                            ...) {
  # return range of sequences for different paths
  # returns NULL if not terminal path is found
  
  seq.list = GetSequenceRanges(grammar, max.depth, startSymb)
  max.seq.len = max(sapply(seq.list, length))
  seq.len.equal = lapply(seq.list, function(x) c(x, rep.int(0, max.seq.len - length(x))))
  max.seq = do.call(pmax, seq.len.equal)
  
  return (as.numeric(max.seq))
}

GetSequenceRanges <- function (grammar, max.depth = length(grammar$def), startSymb = grammar$startSymb) {
  # return range of sequences for different paths
  # returns NULL if not terminal path is found
  
  stopifnot("grammar" %in% class(grammar))
  
  if (IsSymbolTerminal(startSymb)) 
    return(0)
  
  if (max.depth == 0) # this is handled as a "single item"
    return(NULL)
  
  TERMINAL = TRUE
  NON_TERMINAL = FALSE
  
  all_seq = list(NULL)
  while (IsSymbolTerminal(startSymb) == FALSE) {
    sep.symbs = GetFirstNonTerminalandRest(startSymb)
    non.terminal.symb = sep.symbs$nonterminal
    
    # extract the possible sequences
    possible.choices = GetPossibleRuleChoices(non.terminal.symb, grammar)
    seq_list = lapply(1:possible.choices, function(choice.no) {
      chosen.rule = ChosenGrammarRule(non.terminal.symb, choice.no, grammar)
      GetSequenceRanges(grammar, max.depth - 1, startSymb = chosen.rule)
    })
    
    seq_list = unlist(seq_list, recursive=FALSE)
    
    # remove not-terminals
    seq_list = seq_list[!sapply(seq_list, is.null)]
    
    # return NULL if all sequence are non-terminal
    if (length(seq_list) == 0) {
      return(NULL)
    }
    
    # remove the 0's in the list (non-choice elements)
    is.zero = sapply(seq_list, function(x) ifelse(is.numeric(x), x == 0, FALSE))
    seq_list = seq_list[!is.zero]
    if (sum(is.zero) > 0) {
      seq_list[length(seq_list) + 1] = list(NULL)
    }
    
    # append the previously discovered sequence while adding the current choice
    all_seq = lapply(all_seq, function(s) lapply(seq_list, function(s2) c(s, possible.choices, s2)))
    all_seq = unlist(all_seq, recursive=FALSE)
    
    # Do the same thing for the rest
    startSymb = sep.symbs$rest
  }
  return(all_seq)
}

