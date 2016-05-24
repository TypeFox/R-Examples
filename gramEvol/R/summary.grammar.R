summary.grammar <- function(object, ...) {
  ret = list(Start.Symbol = GrammarStartSymbol(object),
             Is.Recursive  = GrammarIsRecursive(object, ...),
             Tree.Depth = GrammarGetDepth(object, ...),
             Maximum.Sequence.Length = GrammarMaxSequenceLen(object, ...),
             Maximum.Rule.Size = GrammarMaxRuleSize(object),
             Maximum.Sequence.Variation = GrammarMaxSequenceRange(object, ...),
             No.of.Unique.Expressions = GrammarNumOfExpressions(object, ...))
  class(ret) <- "summary.grammar"
  
  return(ret)
}

print.summary.grammar <- function(x, ...) {
  cat("Start Symbol:               ", x$Start.Symbol, '\n')
  cat("Is Recursive:               ", x$Is.Recursive, '\n')
  if (x$Is.Recursive)
    cat("Tree Depth:                  Limited to", x$Tree.Depth, '\n')
  else
    cat("Tree Depth:                 ", x$Tree.Depth, '\n')
  cat("Maximum Rule Choices:       ", x$Maximum.Rule.Size, '\n')
  cat("Maximum Sequence Length:    ", x$Maximum.Sequence.Length, '\n')
  cat("Maximum Sequence Variation: ", x$Maximum.Sequence.Variation, '\n')
  cat("No. of Unique Expressions:  ", x$No.of.Unique.Expressions, '\n')
}
