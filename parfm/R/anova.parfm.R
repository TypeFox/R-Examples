################################################################################
#  anova method for class 'parfm'                                              #
################################################################################
#                                                                              #
#  anova.parfm implemets the 'anova' function the objects of class 'parfm'     #
#                                                                              #
#                                                                              #
#   Date: June 5, 2012                                                         #
#   Last modification on: June 5, 2012                                         #
################################################################################

anova.parfm <- function(object,
                        ...) {
  if (length(list(object, ...)) > 1)
    return(anova.parfmlist(list(object, ...)))
  else{
    if (!inherits(object, "parfm"))
      stop(paste("The argument must be a parametric frailty model,",
                 "object of class 'parfm'"))
    
    termlist <- attr(object, "terms")
    
    nmodels <- length(termlist)
    df <- integer(nmodels + 1)
    loglik <- double(nmodels + 1)
    df[1] <- 0
    mycall <- attributes(object)$call
    mycall[["showtime"]] <- FALSE
    mycall[[2]] <- update.formula(mycall[[2]], . ~ 1)
    loglik[1] <- attr(eval(mycall), "loglik")
    for (i in 1:nmodels) {
      df[1 + i] <- df[i] + 1
      mycall[[2]] <- eval(parse(text=paste(
        "update.formula(mycall[[2]], . ~ . +", termlist[i], ")")))
      loglik[1 + i] <- attr(eval(mycall), "loglik")
    }
    table <- data.frame(loglik=loglik, 
                        Chisq=c(NA, 2*diff(loglik)),
                        Df=c(NA, diff(df)))
    table[['Pr(>|Chi|)']] <- 1- pchisq(table$Chisq, table$Df)
    row.names(table) <- c('NULL', termlist)
    title <- paste("Analysis of Deviance Table\n",
                   "Parametric frailty model: response is ",
                   deparse(attr(object, "call")[
                     match("formula", names(attr(object, "call")))][[1]][[2]]),
                   "\nTerms added sequentially (first to last)\n", 
                   sep = "")
    structure(table, heading = title, class = c("anova", "data.frame"))    
  }
}

################################################################################
#  anova method for a list of objects of class 'parfm'                         #
################################################################################
#                                                                              #
#                                                                              #
#   Date: June 5, 2012                                                         #
#   Last modification on: June 5, 2012                                         #
################################################################################

anova.parfmlist <- function(object) {
  if (!is.list(object)) stop("First argument must be a list")
  if (!all(unlist(lapply(object, function(x) inherits(x, 'parfm')))))
    stop("Argument must be a list of 'parfm' models")
  
  nmodels <- length(object)
  if (nmodels == 1) # only one model remains
    return(anova.parfm(object[[1]]))
  
  ### *** Controls *** #########################################################
  responses <- as.character(unlist(lapply(object, function(x) {
    deparse(formula(attr(x, "call"))[[2]])})))
  sameresp <- (responses == responses[1])
  if (!all(sameresp)) {
    object <- object[sameresp]
    warning(paste("Models with response", deparse(responses[!sameresp]), 
                  "removed because response differs from", "model 1"))
  }
  
  whichargs <- names(attr(object[[1]], "call")[
    -match(c("", "formula"), names(attr(object[[1]], "call")))])
  for (arg in whichargs) {
    for (i in 2:nmodels) {
      if (attr(object[[i]], "call")[arg][[1]] != 
        attr(object[[1]], "call")[arg][[1]]) stop(paste(
          "All models must have the same value for argument '",
          arg, "'", sep=""))
    }
  }
  
  ##############################################################################
  
  loglik <- df <- double(nmodels)
  for (i in 1:nmodels) {
    loglik[i] <- attr(object[[i]], "loglik")
    df[i] <- sum(!is.na(object[[i]][, "ESTIMATE"]))
  }
  
  table <- data.frame(loglik, Chisq= c(NA, abs(2*diff(loglik))), 
                      Df= abs(c(NA, diff(df))))
  
  tfun <- function(x) paste(as.character(delete.response(
    terms(formula(attr(x, "call"))))), collapse=' ')
  
  variables <- lapply(object, tfun)
  dimnames(table) <- list(1:nmodels, 
                          c("loglik", "Chisq", "Df"))
  title <- paste("Analysis of Deviance Table\n Parametric frailty model: response is ",
                 responses[1]) 
  topnote <- paste(" Model ", format(1:nmodels), ": ", variables, 
                   sep = "", collapse = "\n")
  table[['P(>|Chi|)']] <- 1-pchisq(table$Chisq, table$Df)
  
  structure(table, heading = c(title, topnote), 
            class = c("anova", "data.frame"))
}
