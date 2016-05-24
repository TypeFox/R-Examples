#----------------------------------------------------------------------------
# localsolver
# Copyright (c) 2014, WLOG Solutions
#----------------------------------------------------------------------------

#' 
#' @demo
#' Google machine reassignment example.
#' Requires license for the full version of localsolver!
#' 

model.text.lsp <- lsp.model.example('extdata/google.txt')

data(google_example_1_data)

lsp <- ls.problem(model.text.lsp)
lsp <- set.params(lsp, 
                  lsTimeLimit=c(150, 50, 40, 30, 30), 
                  indexFromZero=TRUE)

lsp <- add.output.expr(lsp, "x", c(google_example_1_data$nbProcesses, google_example_1_data$nbMachines))
lsp <- add.output.expr(lsp, "totalLoadCost", 1)
lsp <- add.output.expr(lsp, "obj", 1)

ls.solve(lsp, google_example_1_data)

