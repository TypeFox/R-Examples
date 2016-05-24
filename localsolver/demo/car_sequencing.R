#----------------------------------------------------------------------------
# localsolver
# Copyright (c) 2014, WLOG Solutions
#----------------------------------------------------------------------------

#' 
#' @demo
#' The car sequencing example.
#' Requires license for the full version of localsolver!
#' 

model <- "function model() {
  cp[1..nbClasses][1..nbPositions] <- bool();

  // constraints: 
  // for each class c, no more than nbCars[c] assigned to positions
  for [c in 1..nbClasses] 
    constraint sum[p in 1..nbPositions](cp[c][p]) == nbCars[c];

  // constraints: one car assigned to each position p
  for [p in 1..nbPositions] 
    constraint sum[c in 1..nbClasses](cp[c][p]) == 1;

  // expressions: 
  // op[o][p] = 1 if option o appears at position p, and 0 otherwise
  op[o in 1..nbOptions][p in 1..nbPositions] 
    <- or[c in 1..nbClasses : options[c][o]](cp[c][p]);

  // expressions: compute the number of cars in each window
  nbCarsWindows[o in 1..nbOptions][p in 1..nbPositions-ratioDenoms[o]+1] 
    <- sum[k in 1..ratioDenoms[o]](op[o][p+k-1]);

  // expressions: compute the number of violations in each window
  nbViolationsWindows[o in 1..nbOptions][p in 1..nbPositions-ratioDenoms[o]+1]
    <- max(nbCarsWindows[o][p]-ratioNums[o], 0);

  // objective: minimize the sum of violations for all options and all windows
  obj <- sum[o in 1..nbOptions][p in 1..nbPositions-ratioDenoms[o]+1](nbViolationsWindows[o][p]);
  minimize obj;
}"

nbClasses <- 22L
nbPositions <- 100L
nbOptions <- 5L

lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=60, lsTimeBetweenDisplays = 3, lsNbThreads=4)
lsp <- add.output.expr(lsp, expr.text.lsp="cp", dimensions=c(nbClasses, nbPositions))

data <- list(nbPositions = nbPositions,
             nbClasses = nbClasses,
             nbOptions = nbOptions,
             ratioNums = c(1L, 2L, 1L, 2L, 1L),
             ratioDenoms = c(2L, 3L, 3L, 5L, 5L),
             nbCars = c(6L, 10L, 2L, 2L, 8L, 15L, 1L, 5L, 2L, 3L, 2L, 1L, 8L, 3L, 10L, 4L, 4L, 2L, 4L, 6L, 1L, 1L),
             options = matrix(
               data=
                 c(1L, 0L, 0L, 1L, 0L,
                   1L, 1L, 1L, 0L, 0L,
                   1L, 1L, 0L, 0L, 1L,
                   0L, 1L, 1L, 0L, 0L,
                   0L, 0L, 0L, 1L, 0L,
                   0L, 1L, 0L, 0L, 0L,
                   0L, 1L, 1L, 1L, 0L,
                   0L, 0L, 1L, 1L, 0L,
                   1L, 0L, 1L, 1L, 0L, 
                   0L, 0L, 1L, 0L, 0L,
                   1L, 0L, 1L, 0L, 0L, 
                   1L, 1L, 1L, 0L, 1L, 
                   0L, 1L, 0L, 1L, 0L, 
                   1L, 0L, 0L, 1L, 1L, 
                   1L, 0L, 0L, 0L, 0L, 
                   0L, 1L, 0L, 0L, 1L, 
                   0L, 0L, 0L, 0L, 1L, 
                   1L, 0L, 0L, 0L, 1L, 
                   1L, 1L, 0L, 0L, 0L, 
                   1L, 1L, 0L, 1L, 0L, 
                   1L, 0L, 1L, 0L, 1L,
                   1L, 1L, 1L, 1L, 1L
                ), nrow = nbClasses, ncol = nbOptions, byrow = TRUE))


ls.solve(lsp, data)

