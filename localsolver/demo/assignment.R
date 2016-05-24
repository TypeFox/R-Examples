#----------------------------------------------------------------------------
# localsolver
# Copyright (c) 2014, WLOG Solutions
#----------------------------------------------------------------------------

#'
#' @demo  
#' The problem consists in assigning 4 tasks to 4 people. Each person has certain preferences as for the tasks.
#' Each person should get only one task. Each task should be assigned to somebody. The objective is to distribute
#' the tasks so that the sum of tohe preferences of all the people is maximized.
#' 

model <- "function model() {

 //assignments
 assignment[i in 1..4][j in 1..4] <- bool();

 //columns sum
 for [j in 1..4]
  constraint sum[i in 1..4](assignment[i][j]) == 1;

 //rows sum
 for [i in 1..4]
  constraint sum[j in 1..4](assignment[i][j]) == 1;

 //maximize prefernces
 preferencesSum <- sum[i in 1..4][j in 1..4] (preferences[i][j] * assignment[i][j]);
 maximize preferencesSum;

}"

lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=60, lsIterationLimit=250)
lsp <- add.output.expr(lsp, "preferencesSum")
lsp <- add.output.expr(lsp, "assignment", dimensions=c(4,4))

data <- list(preferences=matrix(c(1L,2L,3L,4L,2L,4L,3L,1L,4L,3L,1L,2L,2L,3L,1L,4L), byrow=T,ncol=4)) 
ls.solve(lsp, data)
