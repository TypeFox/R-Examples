cibartControl <- function(n.sim = 20L,
                          n.burn.init = 500L,
                          n.burn.cell = as.integer(n.burn.init / 5L),
                          n.thin = 10L,
                          n.thread = guessNumCores()) {
  for (name in names(formals(cibartControl))) assign(name, as.integer(get(name)))
  
  structure(namedList(n.sim,
                      n.burn.init,
                      n.burn.cell,
                      n.thin,
                      n.thread),
            class = c("cibartControl"))
}

## This assumes that it is 2 frames away from the user, i.e. it has been called by something
## like treatSens.BART or cibart, which has been called (or evaluated such that) in the user's
## environment. If you call it directly, you shouldn't.
evaluateTreatmentModelArgument <- function(arg)
{
  isAnyOf <- function(x, classes) any(sapply(classes, function(class) is(x, class)))
  
  trtCall <- if (!is.null(arg)) arg else formals(cibart)$treatmentModel
  
  if (is.character(trtCall)) trtCall <- parse(text = trtCall)[[1]]
  if (is.call(trtCall)) {
    result <- eval(trtCall, getNamespace("treatSens"), parent.frame(2))
  } else if (is.list(trtCall) && isAnyOf(trtCall, c("probitTreatmentModel", "probitEMTreatmentModel", "bartTreatmentModel"))) {
    result <- trtCall
  } else {
    ## could be that an object was specified, could be just "probit" or "bart" which should be
    ## evaluated in namespace; check the latter first
    result <- tryCatch(eval(call(as.character(trtCall)), getNamespace("treatSens"), parent.frame(2)), error = function(e) e)
    if (is(result, "error"))
      result <- tryCatch(get(as.character(trtCall), parent.frame(2)), error = function(e) e)
    
    if (is(result, "error") || !isAnyOf(result, c("probitTreatmentModel", "probitEMTreatmentModel", "bartTreatmentModel")))
      stop("treatment model of unrecognized type")
  }
  result
}

cibart <- function(Y, Z, X, X.test,
                   zetaY, zetaZ, theta,
                   est.type, treatmentModel = probitEM(),
                   control = cibartControl(), verbose = FALSE)
{
  matchedCall <- match.call()
  
  if (!is(control, "cibartControl")) stop("control must be of class cibartControl; call cibartControl() to create");

  treatmentModel <- evaluateTreatmentModelArgument(matchedCall$treatmentModel)
    
  if (is(treatmentModel, "probitTreatmentModel") && !identical(treatmentModel$family, "flat")) {
    treatmentModel$scale <- rep_len(treatmentModel$scale, ncol(X) + 1L)
  }
  
  .Call("treatSens_fitSensitivityAnalysis",
        Y, Z, X,
        X.test,
        zetaY, zetaZ,
        theta, est.type, treatmentModel,
        control, verbose)
}
