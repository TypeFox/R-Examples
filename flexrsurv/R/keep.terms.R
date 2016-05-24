keep.terms0 <- function (termobj, keepx = NULL, keep.response = FALSE) 
{
    if (is.null(keepx)) 
        termobj
    else {
        drop.terms(termobj, dropx=-keepx, keep.response = keep.response)
    }
}

NPHNLL2NLL.terms <- function(termobj){
  if (!inherits(termobj, "terms")) 
    stop("'termobj' must be a object of class ", dQuote("terms"))
  if(!is.null(attr(termobj, "specials")$NPHNLL)){
    # remove unnecessary options
    argsNLL <- formals("NLL")
    newrhs <- attr(termobj, "term.labels")
    for( i in attr(termobj, "specials")$NPHNLL){
      callNPHNLL <- match.call(NPHNLL, attr(termobj,"variables")[[i+1]])
      indxNLL <- match(names(argsNLL ), names(callNPHNLL))
      callNLL <- callNPHNLL[c(1, indxNLL)]
      callNLL[[1]] <- as.name("NLL")
      newrhs[i-1] <- deparse(callNLL, width.cutoff=500)
    }
    newformula <- reformulate(newrhs, 
                              termobj[[2L]], attr(termobj, "intercept"))
    environment(newformula) <- environment(termobj)
    terms(newformula, specials = names(attr(termobj, "specials")))
  }
  else {
    termobj
  }
}


NPHNLL2NPH.terms <- function(termobj){
  if (!inherits(termobj, "terms")) 
    stop("'termobj' must be a object of class ", dQuote("terms"))
  if(!is.null(attr(termobj, "specials")$NPHNLL)){
    # remove unnecessary options
    argsNPH <- formals("NPH")
    newrhs <- attr(termobj, "term.labels")
    for( i in attr(termobj, "specials")$NPHNLL){
      callNPHNLL <- match.call(NPHNLL, attr(termobj,"variables")[[i+1]])
      indxNPH <- match(names(argsNPH ), names(callNPHNLL))
      callNPH <- callNPHNLL[c(1, indxNPH)]
      callNPH[[1]] <- as.name("NPH")
      newrhs[i-1] <- deparse(callNPH, width.cutoff=500)
    }
    newformula <- reformulate(newrhs, 
                              termobj[[2L]], attr(termobj, "intercept"))
    environment(newformula) <- environment(termobj)
    terms(newformula, specials = names(attr(termobj, "specials")))
  }
  else {
    termobj
  }
}

NPHNLL2Spline.terms <- function(termobj){
  if (!inherits(termobj, "terms")) 
    stop("'termobj' must be a object of class ", dQuote("terms"))
  argsTPS <- formals("TPSplineBasis") 
  argsMS <- formals("MSplineBasis") 
  argsBS <- formals("BSplineBasis")
#  rename args of SplineBasis
  names(argsTPS) <- sub("knots", "Knots_t", names(argsTPS))
  names(argsTPS) <- sub("degree", "Degree_t", names(argsTPS))

  names(argsMS) <- sub("knots", "Knots_t", names(argsMS))
  names(argsMS) <- sub("degree", "Degree_t", names(argsMS))
  names(argsMS) <- sub("keep.duplicates", "keep.duplicates_t", names(argsMS))

  names(argsBS) <- sub("knots", "Knots_t", names(argsBS))
  names(argsBS) <- sub("degree", "Degree_t", names(argsBS))
  names(argsBS) <- sub("keep.duplicates", "keep.duplicates_t", names(argsBS))




  if(!is.null(attr(termobj, "specials")$NPHNLL)){
    listSplineBasis <- list()
    for( i in attr(termobj, "specials")$NPHNLL){
      callNPHNLL <- match.call(NPHNLL, attr(termobj,"variables")[[i+1]])
      if(callNPHNLL$Spline=="tp-spline" ){
        indxS <- match(names(argsTPS ), names(callNPHNLL))
        callSParam <- callNPHNLL[c(1,indxS)]
        names(callSParam) <- sub("Knots_t", "knots", names(callSParam))
        names(callSParam) <- sub("Degree_t", "degree", names(callSParam))
        callSParam[[1]] <- as.name("TPSplineBasis")
        listSplineBasis <- c(listSplineBasis, eval(callSParam))
        callSParam[["type"]] <- "standard"
      }
      if(callNPHNLL$Spline=="tpi-spline" ){
        indxS <- match(names(argsTPS ), names(callNPHNLL))
        callSParam <- callNPHNLL[c(1,indxS)]
        names(callSParam) <- sub("Knots_t", "knots", names(callSParam))
        names(callSParam) <- sub("Degree_t", "degree", names(callSParam))
        callSParam[[1]] <- as.name("TPSplineBasis")
        callSParam[["type"]] <- "standard"
        listSplineBasis <- c(listSplineBasis, eval(callSParam))
      }
      else if(callNPHNLL$Spline=="b-spline" ){
        indxS <- match(names(argsMS ), names(callNPHNLL))
        callSParam <- callNPHNLL[c(1,indxS)]
        names(callSParam) <- sub("Knots_t", "knots", names(callSParam))
        names(callSParam) <- sub("Degree_t", "degree", names(callSParam))
        names(callSParam) <- sub("keep.duplicates_t", "keep.duplicates", names(callSParam))
        callSParam[[1]] <- as.name("MSplineBasis")
        listSplineBasis <- c(listSplineBasis, eval(callSParam))
      }
    }
    listSplineBasis
  }
  else {
    NULL
  }
}


keep.terms <- function (termobj, keepx = NULL, keep.response = FALSE, keep.intercept = FALSE) 
{
    if (is.null(keepx)) {
      if (!inherits(termobj, "terms")) 
        stop("'termobj' must be a object of class ", dQuote("terms"))
      newformula <- reformulate(attr(termobj, "term.labels"), 
                                if (keep.response) 
                                termobj[[2L]]
                                else NULL,
                                if (keep.intercept) 
                                attr(termobj, "intercept")
                                else FALSE)
      environment(newformula) <- environment(termobj)
      terms(newformula, specials = names(attr(termobj, "specials")))
    }
    else {
      if (!inherits(termobj, "terms")) 
        stop("'termobj' must be a object of class ", dQuote("terms"))
      newformula <- reformulate(attr(termobj, "term.labels")[keepx], 
                                if (keep.response) 
                                termobj[[2L]]
                                else NULL,
                                if (keep.intercept) 
                                attr(termobj, "intercept")
                                else FALSE)
      environment(newformula) <- environment(termobj)
      terms(newformula, specials = names(attr(termobj, "specials")))
    }
}


