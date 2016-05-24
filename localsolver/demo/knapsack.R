#----------------------------------------------------------------------------
# localsolver
# Copyright (c) 2014, WLOG Solutions
#----------------------------------------------------------------------------

#' 
#' @demo
#' The problem consists in filling a knapsack with some of available items. Each item has its value and its weight.
#' The overall weight of the items in the knapsack should not exceed a weight bound which is determined for the
#' knapsack. The objective is to maximize the sum of values of the items in the knapsack.
#' 

model <- "function model() {
  x[i in 1..4] <- bool();

  // weight constraint
  knapsackWeight <- sum[i in 1..nbItems](itemWeights[i] * x[i]);
  constraint knapsackWeight <= knapsackBound;
  
  // maximize value
  knapsackValue <- sum[i in 1..nbItems](itemValues[i] * x[i]);
  maximize knapsackValue;
}"

lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=60, lsIterationLimit=250)
lsp <- add.output.expr(lsp, "x", 4)
lsp <- add.output.expr(lsp, "knapsackWeight")
lsp <- add.output.expr(lsp, "knapsackValue")

data <- list(nbItems=4L, itemWeights=c(1L,2L,3L,4L), itemValues=c(5,6,7,8), knapsackBound = 9L)
ls.solve(lsp, data)

