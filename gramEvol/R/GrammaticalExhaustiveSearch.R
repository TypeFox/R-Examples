GrammaticalExhaustiveSearch <- function(grammar, evalFunc,
                                        max.depth = GrammarGetDepth(grammar),
                                        startSymb = GrammarStartSymbol(grammar),
                                        max.len = GrammarMaxSequenceLen(grammar, max.depth, startSymb),
                                        wrappings = 3,
                                        terminationCost = NA,
                                        monitorFunc = NULL) {
  
  best.expr = NULL
  best.score = Inf
  best.seq = NULL
  
  iterations = 0
  genome = NULL
  while (TRUE) {
    genome = GrammarGetNextSequence(grammar, genome, startSymb, max.depth, max.len)
    
    if (is.GrammarOverflow(genome))
      break
    
    iterations = iterations + 1
    
    expr = as.expression(GrammarMap(genome, grammar, wrappings))
    
    score = evalFunc(expr)
    
    if (score < best.score) {
      best.score = score
      best.expr = expr
      best.seq = genome
      if (!is.na(terminationCost)) {
        if (score <= terminationCost) {
          break        
        }
      }
    }
    
    # call the monitor function
    if (!is.null(monitorFunc)) {

      res <- list(currentSequence = genome,
                  currentExpression = expr,
                  currentCost = score,
                  bestSequence = best.seq,
                  bestExpression = best.expr,
                  bestCost = best.score,
                  numExpr = iterations)

      class(res) <- "GESearch"
      monitorFunc(res)
    }
  }
  
  res <- list(bestSequence = best.seq,
              bestExpression = best.expr,
              bestCost = best.score,
              numExpr = iterations)
  class(res) <- "GESearch"

  return(res)
}