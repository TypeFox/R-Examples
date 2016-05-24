fix.flexrsurv.formula <- function(formula, data,
                                  na.action=options("na.action"),
                                  contrasts=NULL,
                                  debug=FALSE,
                                  model=c("additive","multiplicative"),
                                  Spline=c("b-spline", "tp-spline", "tpi-spline"),
                                  method =c("glm","MLE")){
# fix the formula for flexrsurv
#

  model  <- match.arg(model)          # type of model used (additive for remontet's model - multiplicative for mahboubi's model)
  method <- match.arg(method)         # method of relative surv. curve estimation
  Spline <- match.arg(Spline)         # choice of spline basis
  
  special <- c("NPH","NLL", "NPHNLL") 
  Terms <- if (missing(data)){
    terms(formula, special)
  } else {
    terms(formula, special, data = data)
  }
  
  
  NamesLinVars <- all_LIN_vars(Terms)
  
  NamesNLLVars <- all_specials_vars(Terms,
                                    specials="NLL",
                                    unique = TRUE,
                                    order="formula")
  NamesNPHVars<- all_specials_vars( Terms,
                                    specials="NPH",
                                    unique = TRUE,
                                    order="formula")
  
  NamesNPHNLLVars<- all_specials_vars( Terms,
                                       specials="NPHNLL",
                                       unique = TRUE,
                                       order="formula")
  
  
  
  modified <- 0
  newtermlabels <- labels(Terms)

  # when a variable is both (either LL or NLL) and NPH, then Intercept.t must be FALSE in NPH()


  if (sum(NamesNPHVars %in% c(NamesLinVars, NamesNLLVars))){
    # some LIN or NLL variables are also NPH
    # force Intercept.t=FALSE in NPH()
    LNPHVars <- NamesNPHVars[NamesNPHVars %in% c(NamesLinVars, NamesNLLVars)]
    
    for (i in attr(Terms, "specials")[["NPH"]]){
      thecall <-  match.call(NPH, attr(Terms,"variables")[[i+1]])
        # the variable name is the second argument of the special function
      if (as.character(thecall[["x"]]) %in% LNPHVars) {
        if ( (length(thecall[["Intercept.t"]]) == 0 ) || ( thecall[["Intercept.t"]] == TRUE ) ){
          thecall[["Intercept.t"]] <- FALSE
          msg1 <- gettextf("both {linear or NPH()} and NPH() effect in the formula for variable %s.\n",
                           sQuote(as.character(thecall[["x"]])),
                           domaine=NA)
          msg2 <- gettextf("Force Intercept.t=FALSE in NPH(%s).", c(thecall[[2]]), domaine=NA)
          warning(paste(msg1, msg2, sep=""), call. = FALSE, immediate. = debug)
          }
          modified <- modified + 1
            # change Intercept.t to FALSE if not set to false
          thecall[["Intercept.t"]] <- FALSE
          indxterm <- variable2term(i, Terms)
          print(indxterm)
          charcall<-deparse(thecall, 500)
          oldtermlabel <- newtermlabels[indxterm[1]]
          newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
      }
    }
  }

  if(modified > 0){
    formula <- reformulate(newtermlabels,
                           response = if (attr(Terms, "response")){ 
                             Terms[[2L]]
                           }
                           else NULL,
                           intercept = attr(Terms, "intercept"))
  }
  
  
  Terms <- if (missing(data)){
    terms(formula, special)
  } else {
    terms(formula, special, data = data)
  }
  
  modified <- 0
  newtermlabels <- labels(Terms)
  # get the names of the end of intervals
  thesurvcall <-  match.call(Surv, attr(Terms,"variables")[[2]])
  # print(thesurvcall)
  survtimename <- if(length(thesurvcall[["event"]])>0){
    as.character(thesurvcall[["time2"]])
  } else {
    as.character(thesurvcall[["time"]])
  }
  
  # force type="right" if time2 is given
  #  if(length(thesurvcall[["event"]])>0){
  #    if(length(thesurvcall[["type"]])==0){
  #      thesurvcall[["type"]] <- "right"
  #      formula <- reformulate(labels(Terms),response = thesurvcall,intercept = attr(Terms, "intercept"))
  #    }
  #  }
  
  
  # check that NPH() time variable is the same as Surv() time
  if(length(NamesNPHVars) >0){
    for (i in attr(Terms, "specials")[["NPH"]]){
      thecall <-  match.call(NPH, attr(Terms,"variables")[[i+1]])
      timename <- as.character(thecall[["timevar"]])
      if (timename != survtimename) {
        stop("The time variable in NPH effect must be the same as in the survival object.\n")
      }
    }
  }
  
  
  # check that NPHNLL() time variable is the same as Surv() time
  if(length(NamesNPHNLLVars) >0){
    for (i in attr(Terms, "specials")[["NPHNLL"]]){
        thecall <-  match.call(NPHNLL, attr(Terms,"variables")[[i+1]])
        timename <- as.character(thecall[["timevar"]])     
        if (timename != survtimename) {
          stop("the time variable in NPHNLL effect must be the same as in the survival object. \n")
        }
      }
  }
  # add model type (additive ormultiplicative)
#  if (length(grep("model", formula, value=TRUE)) & method == "glm") {
#    warning(paste("The model option in NPHNLL() effects is overwritten by the flexrsurv(model=", model, ") arguments.", sep=""),
#            call. = FALSE#, immediate. = debug
#            )
#  }
  if(length(NamesNPHNLLVars) >0){
    for (i in attr(Terms, "specials")[["NPHNLL"]]){
      thecall <-  match.call(NPHNLL, attr(Terms,"variables")[[i+1]])
      if ( length(thecall[["model"]])==0 ){
        modified <- modified + 1
        thecall[["model"]] <- model
        indxterm <- variable2term(i, Terms)
        charcall<-deparse(thecall, 500)
        oldtermlabel <- newtermlabels[indxterm[1]]
        newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
      } else if ((thecall[["model"]] != model) & (method == "glm")){
        warning(gettextf("The 'model' option in NPHNLL() effects is overwritten by model=%s.\n", dQuote(model), domaine=NA),
                call. = FALSE)
        modified <- modified + 1
        thecall[["model"]] <- model
        indxterm <- variable2term(i, Terms)
        charcall<-deparse(thecall, 500)
        oldtermlabel <- newtermlabels[indxterm[1]]
        newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
      }
    }
  }
  
  
  if(modified > 0){
    formula <- reformulate(newtermlabels,
                           response = if (attr(Terms, "response")){ 
                             Terms[[2L]]
                           }
                           else NULL,
                           intercept = attr(Terms, "intercept"))
  }
  
  
  Terms <- if (missing(data)){
    terms(formula, special)
  } else {
    terms(formula, special, data = data)
  }
  
  modified <- 0
  newtermlabels <- labels(Terms)
  
  if (length(grep("Spline", formula, value=TRUE))) {
    warning(gettextf("The 'Spline' option in NLL(), NPH() and NPHNLL() effects is overwritten by Spline=%s.\n",
                  dQuote(Spline), domaine=NA),
            call. = FALSE#, immediate. = debug
            )
  }
  
  
  # set type of splines
  # evaluate degrees à knots & boundaries
  if(length(NamesNLLVars) >0){
    for (i in attr(Terms, "specials")[["NLL"]]){
      thecall <-  match.call(NLL, attr(Terms,"variables")[[i+1]])
      if ( length(thecall[["Spline"]])==0 ){
        modified <- modified + 1
        thecall[["Spline"]] <- Spline
      }
      if ( length(thecall[["Degrees"]])!=0 ){
        thecall[["Degrees"]] <- eval(as.expression(thecall[["Degrees"]]))
        modified <- modified + 1
      }
      if ( length(thecall[["Knots"]])!=0 ){
        thecall[["Knots"]] <- eval(as.expression(thecall[["Knots"]]))
        modified <- modified + 1
      }
      if ( length(thecall[["Boundary.knots"]])!=0 ){
        thecall[["Boundary.knots"]] <- eval(as.expression(thecall[["Boundary.knots"]]))
        modified <- modified + 1
      }
      indxterm <- variable2term(i, Terms)
      charcall<-deparse(thecall, 500)
      oldtermlabel <- newtermlabels[indxterm[1]]
      newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
    }
  }
  if(length(NamesNPHVars) >0){
    for (i in attr(Terms, "specials")[["NPH"]]){
      thecall <-  match.call(NPH, attr(Terms,"variables")[[i+1]])
      if ( length(thecall[["Spline"]])==0 ){
        modified <- modified + 1
        thecall[["Spline"]] <- Spline
      }
      if ( length(thecall[["Degrees"]])!=0 ){
        thecall[["Degrees.t"]] <- eval(as.expression(thecall[["Degrees.t"]]))
        modified <- modified + 1
      }
      if ( length(thecall[["Knots"]])!=0 ){
        thecall[["Knots.t"]] <- eval(as.expression(thecall[["Knots.t"]]))
        modified <- modified + 1
      }
      if ( length(thecall[["Boundary.knots.t"]])!=0 ){
        thecall[["Boundary.knots.t"]] <- eval(as.expression(thecall[["Boundary.knots.t"]]))
        modified <- modified + 1
      }
      indxterm <- variable2term(i, Terms)
      charcall<-deparse(thecall, 500)
      oldtermlabel <- newtermlabels[indxterm[1]]
      newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
    }
  }
  
  if(length(NamesNPHNLLVars) >0){
    for (i in attr(Terms, "specials")[["NPHNLL"]]){
      thecall <-  match.call(NPHNLL, attr(Terms,"variables")[[i+1]])
      if ( length(thecall[["Spline"]])==0 ){
        modified <- modified + 1
        thecall[["Spline"]] <- Spline
      }
      if ( length(thecall[["Degrees"]])!=0 ){
        thecall[["Degrees"]] <- eval(as.expression(thecall[["Degrees"]]))
        modified <- modified + 1
      }
      if ( length(thecall[["Knots"]])!=0 ){
        thecall[["Knots"]] <- eval(as.expression(thecall[["Knots"]]))
        modified <- modified + 1
      }
      if ( length(thecall[["Boundary.knots"]])!=0 ){
        thecall[["Boundary.knots"]] <- eval(as.expression(thecall[["Boundary.knots"]]))
        modified <- modified + 1
      }
      if ( length(thecall[["Degrees"]])!=0 ){
        thecall[["Degrees.t"]] <- eval(as.expression(thecall[["Degrees.t"]]))
        modified <- modified + 1
      }
      if ( length(thecall[["Knots"]])!=0 ){
        thecall[["Knots.t"]] <- eval(as.expression(thecall[["Knots.t"]]))
        modified <- modified + 1
      }
      if ( length(thecall[["Boundary.knots.t"]])!=0 ){
        thecall[["Boundary.knots.t"]] <- eval(as.expression(thecall[["Boundary.knots.t"]]))
        modified <- modified + 1
      }
      indxterm <- variable2term(i, Terms)
      charcall<-deparse(thecall, 500)
      oldtermlabel <- newtermlabels[indxterm[1]]
      newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
    }
  }
  if(modified > 0){
    #.#print(Terms[[2L]])
    formula <- reformulate(newtermlabels,
                           response = if (attr(Terms, "response")){ 
                             Terms[[2L]]
                           }
                           else NULL,
                           intercept = attr(Terms, "intercept"))
  }


return(formula)
}
