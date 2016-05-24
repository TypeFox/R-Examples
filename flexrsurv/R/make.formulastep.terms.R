make.formulastepNLL.terms <- function(terms, data,
                                      response=".fail",
                                      baseline= "NLL(.t, Spline = \"b-spline\", Knots = NULL, Degree = 2, Log = FALSE, Intercept = TRUE)",
                                      tik="tik",
                                      ...){
# make formula for glm() at step NLL


  NamesNPHNLLVars<- all_specials_vars( terms,
                                      specials="NPHNLL",
                                      unique = TRUE,
                                      order="formula")
  

  modified <- 0
  newtermlabels <- labels(terms)
  
      # add offset(log(tik))
  offset <- paste("offset(log(", tik, "))", sep="")
  
      # add arg BETAt = "betaT"x" in NPHNLL() call
  # force intercept.t = FALSE
  if(length(NamesNPHNLLVars) >0){
    for (i in attr(terms, "specials")["NPHNLL"]){
#      for (k in 1:length(i)){
        thecall <-  match.call(NPHNLL, attr(terms,"variables")[[i+1]])
        namebetaTx <- paste("betaT", thecall[["x"]], sep="")
        modified <- modified + 1
        thecall[[1]] <- as.name("NLLbeta")
        thecall[["y"]] <- as.name(namebetaTx)
        indxterm <- variable2term(i, terms)
        charcall<-deparse(thecall, 500)
        oldtermlabel <- newtermlabels[indxterm]
        newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
        if(thecall[["Spline"]]=="b-spline"){
          # add offset(alpha(x)*b1
          minX <- eval(as.expression(thecall[["Boundary.knots"]]))[1]
#          offset <- c(offset, paste("offset(", namebetaTx, " * ", minX, ")", sep=""))
        }
      }
#    }
  }
  
  if(modified > 0){
    formula <- reformulate(c(baseline, newtermlabels, offset), 
                           response = response,
                           intercept = FALSE)
  }
  
return(formula)
  
}




make.formulastepNPH.terms <- function(terms, data,
                                       response=".fail",
                                       baseline= "NLL(.t, Spline = \"b-spline\", Knots = NULL, Degree = 2, Log = FALSE, Intercept = TRUE)",
                                       tik="tik",
                                       ...){
# make formula for glm() at step NPH


  
  NamesNPHNLLVars<- all_specials_vars( terms,
                                      specials="NPHNLL",
                                      unique = TRUE,
                                      order="formula")
  

  modified <- 0
  newtermlabels <- labels(terms)

      # add offset(log(tik))
  offset <- paste("offset(log(", tik, "))", sep="")
  
      # change change arg x to alpha"x" in NPHNLL() call
  # force intercept.t = FALSE
  if(length(NamesNPHNLLVars) >0){
    for (i in attr(terms, "specials")["NPHNLL"]){
        thecall <-  match.call(NPHNLL, attr(terms,"variables")[[i+1]])
        namealphax <- paste("alpha", thecall[["x"]], sep="")
        modified <- modified + 1
        thecall[[1]] <- as.name("NPH")
        thecall[["x"]] <- as.name(namealphax)
        thecall[["Intercept.t"]] <- FALSE
        indxterm <- variable2term(i, terms)
        charcall<-deparse(thecall, 500)
        oldtermlabel <- newtermlabels[indxterm]
        newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
        # add offset(alpha(x)*b1
        offset <- c(offset, paste("offset(", namealphax, "b1", ")", sep=""))
      }
  }
  if(modified > 0){
    formula <- reformulate(c(baseline, newtermlabels, offset), 
                           response = response,
                           intercept = FALSE)
  }
  
return(formula)
  
}
