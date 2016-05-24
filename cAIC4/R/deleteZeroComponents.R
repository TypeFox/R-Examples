deleteZeroComponents <-
function(m) {
  # A function that deletes all random effects terms if corresponding variance
  # parameter is estimated to zero.
  #
  # Args: 
  #   m     = Object of class lmerMod. Obtained by lmer()
  #
  # Returns:
  #   m/newMod = A model without zero estimated variance component
  #
  theta      <- getME(m, "theta")
  thetazero  <- which(theta == 0)
  
  if (length(thetazero) == 0) {  # every thing's fine
    return(m)
  }
    
  if (length(theta) == length(thetazero)) {  #  only lm left
    warning("Model has no random effects variance components larger than zero.")
    return(lm(nobars(formula(m)), model.frame(m)))
  }

  varBlockMatrices <- getME(m, "ST")
  cnms <- m@cnms
  
  if(exists("gamm4", m@optinfo)) {  # for gamm4 what to exclude from the model
    for(i in 1:length(varBlockMatrices)){
      if(any(diag(varBlockMatrices[[i]]) == 0)) {
         cat("The term", cnms[[i]][which(diag(varBlockMatrices[[i]]) == 0)], 
          "has zero variance components. \n")
      }
    }
    stop("After removing the terms with zero variance components and refitting 
          the model cAIC can be called again.", call. = FALSE)
  }

  if(is.null(m@optinfo$conv$lme4$code)) {
    for(i in 1:length(varBlockMatrices)){
      cnms[[i]] <- cnms[[i]][which(diag(varBlockMatrices[[i]]) != 0)]
    }
  } else {    # in case of convergence failures
    nc  <- vapply(cnms, length, 1L)
    thl <- split(theta, rep.int(seq_along(nc), (nc * (nc + 1))/2))
    for (i in 1:length(nc)) {
      ranVars   <- thl[[i]][1:nc[i]]
      cnms[[i]] <- cnms[[i]][which(ranVars != 0)] 
    }    
  }
  
  reFormula  <- cnms2formula(cnms)
  if(nobars(formula(m)) == formula(m)[[2]]) {  # if there are no fixed effects 
    rhs      <- reFormula
  } else {
    rhs      <- c(attr(terms(nobars(formula(m))), "term.labels"), reFormula)
  }
  lhs        <- formula(m)[[2]]  # left hand side of the formula
  newFormula <- reformulate(rhs, lhs)  # merge both sides           
  newMod     <- update(m, formula. = newFormula, evaluate = TRUE)
  
  return(deleteZeroComponents(newMod))
}
