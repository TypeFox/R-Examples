library("gramEvol")

ruleDef <- list(
  expr     = grule(coeff * expr, coeff * var, var, var + var),
  var      = grule(A, B),
  coeff    = grule(c1, c2, c3))

# Creating the grammar object
grammar <- CreateGrammar(ruleDef) 

max.range <- GrammarMaxSequenceRange(grammar)
stopifnot(all(max.range == c(4,3,4,3,4,3,2)))

num.expr <- GrammarNumOfExpressions(grammar)
stopifnot(num.expr == 156)

# first sequence
genome <- GrammarGetNextSequence(grammar)
stopifnot(all(genome == c(0,0,0,0,1,0,0)))

genome <- GrammarGetNextSequence(grammar, genome)
stopifnot(all(genome == c(0,0,0,0,1,0,1)))

# start in a loop
genome <- NULL
cnt <- 0
for (i in 1:num.expr) {
  genome <- GrammarGetNextSequence(grammar, genome)
  stopifnot(!is.GrammarOverflow(genome))
  stopifnot(genome < max.range[1:length(genome)])
  cnt <- cnt + 1
}
stopifnot(all(genome == c(3,1,1)))
stopifnot(cnt == num.expr)

# next one should be overflow
genome <- GrammarGetNextSequence(grammar, genome)
stopifnot(is.GrammarOverflow(genome))
