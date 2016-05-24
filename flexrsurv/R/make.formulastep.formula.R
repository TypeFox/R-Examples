make.formulastepNLL.formula <- function(formula, data,
                                      response=".fail",
                                      baseline= "NLL(.t, Spline = \"b-spline\", Knots = NULL, Degree = 2, Log = FALSE, Intercept = TRUE)",
                                      tik="tik",
                                      ...){
# make formula for glm() at step NLL

  special <- c("NPH","NLL", "NPHNLL") 
  Terms <- if (missing(data)){
    terms(formula, special, keep.order = TRUE)
  } else {
    terms(formula, special, data = data, keep.order = TRUE)
  }


  
  NamesNPHNLLVars<- all_specials_vars( Terms,
                                      specials="NPHNLL",
                                      unique = TRUE,
                                      order="formula")
  

  modified <- 0
  newtermlabels <- labels(Terms)
  
      # add offset(log(tik))
  offset <- paste("offset(log(", tik, "))", sep="")
  
      # add arg BETAt = "betaT"x" in NPHNLL() call
  # force intercept.t = FALSE
  if(length(NamesNPHNLLVars) >0){
    for (i in attr(Terms, "specials")["NPHNLL"]){
      for (k in 1:length(i)){
        thecall <-  match.call(NPHNLL, attr(Terms,"variables")[[i[k]+1]])
        namebetaTx <- paste("betaT", thecall[["x"]], sep="")
        modified <- modified + 1
        thecall[[1]] <- as.name("NLLbeta")
        thecall[["y"]] <- as.name(namebetaTx)
        indxterm <- variable2term(i, Terms)
        charcall<-deparse(thecall, 500)
        oldtermlabel <- newtermlabels[indxterm[k]]
        newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
        if(thecall[["Spline"]]=="b-spline"){
          # add offset(alpha(x)*b1
          minX <- eval(as.expression(thecall[["Boundary.knots"]]))[1]
          offset <- c(offset, paste("offset(", namebetaTx, " * ", minX, ")", sep=""))
        }
      }
    }
  }
  
  if(modified > 0){
    formula <- reformulate(c(baseline, newtermlabels, offset), 
                           response = response,
                           intercept = FALSE)
  }
  
  
  
return(formula)
  
}


################################################################################

make.formulastepNPH.formula <- function(formula, data,
                                       response=".fail",
                                       baseline= "NLL(.t, Spline = \"b-spline\", Knots = NULL, Degree = 2, Log = FALSE, Intercept = TRUE)",
                                       tik="tik",
                                       ...){
# make formula for glm() at step NLL

  special <- c("NPH","NLL", "NPHNLL") 
  Terms <- if (missing(data)){
    terms(formula, special, keep.order = TRUE)
  } else {
    terms(formula, special, data = data, keep.order = TRUE)
  }


  
  NamesNPHNLLVars<- all_specials_vars( Terms,
                                      specials="NPHNLL",
                                      unique = TRUE,
                                      order="formula")
  

  modified <- 0
  newtermlabels <- labels(Terms)

      # add offset(log(tik))
  offset <- paste("offset(log(", tik, "))", sep="")
  
      # change change arg x to alpha"x" in NPHNLL() call
  # force intercept.t = FALSE
  if(length(NamesNPHNLLVars) >0){
    for (i in attr(Terms, "specials")["NPHNLL"]){
      for (k in 1:length(i)){
        thecall <-  match.call(NPHNLL, attr(Terms,"variables")[[i[k]+1]])
        namealphax <- paste("alpha", thecall[["x"]], sep="")
        modified <- modified + 1
        thecall[[1]] <- as.name("NPH")
        thecall[["x"]] <- as.name(namealphax)
        thecall[["Intercept.t"]] <- FALSE
        indxterm <- variable2term(i, Terms)
        charcall<-deparse(thecall, 500)
        oldtermlabel <- newtermlabels[indxterm[k]]
        newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
        # add offset(alpha(x)*b1
        offset <- c(offset, paste("offset(", namealphax, "b1", ")", sep=""))
      }
    }
  }
  
  if(modified > 0){
    formula <- reformulate(c(baseline, newtermlabels, offset), 
                           response = response,
                           intercept = FALSE)
  }
  
  
  
return(formula)
  
}


######################################################################

NPHNLL2NLL.formula <- function(formula, data,
                                      response=".fail",
                                      ...){
# make formula in which NPHNLL is replace by NLL

  special <- c("NPH","NLL", "NPHNLL") 
  Terms <- if (missing(data)){
    terms(formula, special, keep.order = TRUE)
  } else {
    terms(formula, special, data = data, keep.order = TRUE)
  }


  
  NamesNPHNLLVars<- all_specials_vars( Terms,
                                      specials="NPHNLL",
                                      unique = TRUE,
                                      order="formula")
  

  modified <- 0
  newtermlabels <- labels(Terms)
  
      # add arg BETAt = "betaT"x" in NPHNLL() call
  # force intercept.t = FALSE
  if(length(NamesNPHNLLVars) >0){
    for (i in attr(Terms, "specials")["NPHNLL"]){
      thecall <-  match.call(NPHNLL, attr(Terms,"variables")[[i+1]])
      modified <- modified + 1
      thecall[[1]] <- as.name("NLL")
      thecall[["timevar"]] <- NULL
      if(!is.null(thecall[["Degree.t"]])){
        thecall[["Degree.t"]] <- NULL
      }
      if(!is.null(thecall[["Knots.t"]])){
        thecall[["Knots.t"]] <- NULL
      }
      if(!is.null(thecall[["Intercept.t"]])){
        thecall[["Intercept.t"]] <- NULL
      }
      if(!is.null(thecall[["Boundary.knots.t"]])){
        thecall[["Boundary.knots.t"]] <- NULL
      }
      if(!is.null(thecall[["model"]])){
        thecall[["model"]] <- NULL
      }
      indxterm <- variable2term(i, Terms)
      charcall<-deparse(thecall, 500)
      oldtermlabel <- newtermlabels[indxterm]
      newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
    }
  }
  if(modified > 0){
    formula <- reformulate(newtermlabels, 
                           response = response,
                           intercept = FALSE)
  }
  
return(formula)
  
}


################################################################################

NPHNLL2NPHalpha.formula <- function(formula, data,
                                       response=".fail",
                                       ...){
# make formula in which NPHNLL is replace by NPHalpha

  special <- c("NPH","NLL", "NPHNLL") 
  Terms <- if (missing(data)){
    terms(formula, special, keep.order = TRUE)
  } else {
    terms(formula, special, data = data, keep.order = TRUE)
  }


  
  NamesNPHNLLVars<- all_specials_vars( Terms,
                                      specials="NPHNLL",
                                      unique = TRUE,
                                      order="formula")
  

  modified <- 0
  newtermlabels <- labels(Terms)
  
      # change change arg x to alpha"x" in NPHNLL() call
  # force intercept.t = FALSE
  if(length(NamesNPHNLLVars) >0){
    for (i in attr(Terms, "specials")["NPHNLL"]){
      thecall <-  match.call(NPHNLL, attr(Terms,"variables")[[i+1]])
      modified <- modified + 1
      thecall[[1]] <- as.name("NPHalpha")
      if(!is.null(thecall[["Degree"]])){
        thecall[["Degree"]] <- NULL
      }
      if(!is.null(thecall[["Knots"]])){
        thecall[["Knots"]] <- NULL
      }
      if(!is.null(thecall[["Intercept"]])){
        thecall[["Intercept"]] <- NULL
      }
      if(!is.null(thecall[["Boundary.knots"]])){
        thecall[["Boundary.knots"]] <- NULL
      }
      indxterm <- variable2term(i, Terms)
      charcall<-deparse(thecall, 500)
      oldtermlabel <- newtermlabels[indxterm]
      newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
    }
  }
  if(modified > 0){
    formula <- reformulate(newtermlabels, 
                           response = response,
                           intercept = FALSE)
  }
  
  
  
return(formula)
  
}
