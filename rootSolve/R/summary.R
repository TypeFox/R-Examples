### ============================================================================
### Summaries of ode variables
### ============================================================================
summary.rootSolve <- function(object, ...){

  att <- attributes(object)
  nspec  <- att$nspec          # for models solved with ode.1D, ode.2D
  dimens <- att$dimens
  svar   <- nspec*prod(dimens)   # number of state variables
  if (is.null(nspec)) nspec <- 1
  if (is.null(dimens)) dimens <- length(object$y)/nspec
  
  # variable names: information for state and ordinary variables is different
  if (is.null(att$ynames))
    if (is.null(dimens))
      varnames <- colnames(object)[1:nspec]
    else
      varnames <- 1:nspec
  else
    varnames <- att$ynames   # this gives one name for multi-dimensional var.

  varnames <- c(varnames, names(object[-1]))
  # summaries for all variables
  Summ <- NULL

  # first state variables
  for (i in 1:nspec) {
    out  <- object[[1]] [((i-1)*prod(dimens)+1):(i*prod(dimens))]
    Summ <- rbind(Summ, c(summary(out, ...), N = length(out), sd = sd(out)))
  }
  # then oridinary variables
  il <- length(object[-1])
  if (il > 0) 
   for (i in 1 : il) {
    out <- as.vector(object[[i+1]])
    Summ <- rbind(Summ, c(summary(out, ...), N = length(out), sd = sd(out)))
  }  

  rownames(Summ) <- varnames  # rownames or an extra column?
  data.frame(t(Summ))         # like this or not transposed?
}

