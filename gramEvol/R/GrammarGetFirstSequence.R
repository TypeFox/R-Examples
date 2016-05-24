GrammarGetFirstSequence <- function (grammar, 
                                     seqStart = NULL, 
                                     startSymb = GrammarStartSymbol(grammar),
                                     max.depth = GrammarGetDepth(grammar),
                                     max.len = GrammarMaxSequenceLen(grammar, max.depth, startSymb)) {
  
  GrammarGetFirstSequence.Recursive(grammar, 
                                    seqStart = seqStart, 
                                    startSymb = startSymb, 
                                    max.depth = max.depth, 
                                    max.len = max.len, 
                                    depth.map = NULL)
}


UpdateDepthMap <- function(depth.map, newSymb, oldSymb) {
  # define the new depth based on how many symbols were replaced
  numOldNonTerms = GetNumNonTerminals(oldSymb)
  numNewNonTerms = GetNumNonTerminals(newSymb)
  numAddedNonTerms = numNewNonTerms - numOldNonTerms
  lenDepth = length(depth.map)
  
  # coming out of len, as one symbol is replaced
  if (lenDepth > 0) {
    depth.map[lenDepth] = depth.map[lenDepth] - 1
  }
  
  # adding to the depth if symbols are added
  if (numAddedNonTerms >= 0) {
    depth.map = c(depth.map, numAddedNonTerms + 1)
  }
  
  # clean all depths that are zero or less
  while(length(depth.map) > 0) {
    if (tail(depth.map, 1) <= 0) {
      depth.map = depth.map[-length(depth.map)]
    } else {
      break
    }
  }
  
  return (depth.map)
}

GrammarGetFirstSequence.Recursive <- function (grammar, 
                                     seqStart = NULL, 
                                     startSymb = GrammarStartSymbol(grammar),
                                     max.depth = GrammarGetDepth(grammar),
                                     max.len = GrammarMaxSequenceLen(grammar, max.depth, startSymb),
                                     depth.map = NULL) {
  
  if (IsSymbolTerminal(startSymb)) {
    return (NULL)
  }
  if (max.len == 0){
    return (NewGrammarOverFlow())
  }
  if  (max.depth - length(depth.map) <= 0) {
    return (NewGrammarOverFlow())
  }
  if (GetNumNonTerminals(startSymb) > max.len)  {
    return (NewGrammarOverFlow())
  }
  
  # find the possible choices
  possible.choices = GetPossibleRuleChoicesFirstSymbol(startSymb, grammar)
  
  # extract the first element
  if (length(seqStart) == 0) {
    seqStart = 0
  }
  current.codon = seqStart[1]
  restOfSequence = seqStart[-1]
  
  # check for overflow
  if (current.codon >= possible.choices) {
    return(NewGrammarOverFlow())
  }
  
  for (i in current.codon:(possible.choices-1)) {
    # apply the rule
    newSymb = ApplyGrammarRule(i, startSymb, grammar)

    # get the new depth map
    new.depth.map = UpdateDepthMap(depth.map, newSymb, oldSymb = startSymb)
    
    # recursively apply for the next
    new_seq = GrammarGetFirstSequence.Recursive(grammar, seqStart = restOfSequence, startSymb = newSymb, 
                                      max.len = max.len - 1, 
                                      max.depth = max.depth,
                                      depth.map = new.depth.map)
    
    if (!is.GrammarOverflow(new_seq)) {
      new_seq = c(i, new_seq)
      break
    }
    
    # if we are past here, rest of Sequence is useless
    restOfSequence = NULL
  }
  
  return (new_seq)
}
