## $Id: valid.R 5596 2014-07-16 22:17:15Z bpikouni $ UTC

## Functions for verifying input arguments to functions

validAlpha <- function(x) {
  ##
  ## PURPOSE: Alpha represents the probability of a false positive in
  ## statistical hypothesis testing (Type I error) and needs to
  ## follow the condition below.
  x <- as.numeric(x)
  if(!(is.numeric(x) && (length(x)==1) && (x > 0) && (x < 1))) {
    stop(cgMessage("alpha needs to be a strictly numeric single",
                   "value that is greater than 0",
                   "and less than 1. Traditional practice uses the",
                   "default setting of 0.05."
                   ))
  }
  else return(TRUE)
}

validPower <- function(x) {
  ##
  ## PURPOSE: Power represents the probability of a true positive in
  ## statistical hypothesis testing and needs to
  ## follow the condition below.
  x <- as.numeric(x)
  if(!(is.numeric(x) && (length(x)==1) && (x > 0) && (x < 1))) {
    stop(cgMessage("Power needs to be a strictly numeric single",
                   "value that is greater than 0",
                   "and less than 1. Traditional practice uses the",
                   "default setting of 0.80."
                   ))
  }
  else return(TRUE)
}

validBoolean <- function(arg) {
  ##
  ## PURPOSE: Simple check that a valid boolean value was specified
  ##
  argname <- deparse(substitute(arg))

  if(!is.logical(arg)) {
    stop(cgMessage("The argument ", argname, " needs to be",
                   "set to boolean TRUE or FALSE.",
                   "The current value seems to be ", arg, "."
                   ))
  }
  else return(TRUE)
}


validCharacter <- function(arg) {
  ##
  ## PURPOSE: Simple check that a valid boolean value was specified
  ##
  argname <- deparse(substitute(arg))

  if(!is.character(arg)) {
    stop(cgMessage("The argument ", argname, " needs to be",
                   "of type character.",
                   "The current value seems to be ", arg, "."
                   ))
  }
  else return(TRUE)
}


validNumeric <- function(arg, positive=FALSE, integer=FALSE) {
  ##
  ## PURPOSE: Simple check that a numeric value was specified
  ##
  if(missing(arg)) {
      stop(cgMessage("The argument ",
                     deparse(substitute(arg)),
                     " is missing and need to be specified."))
  }
  argname <- deparse(substitute(arg))
  validBoolean(positive)
  validBoolean(integer)

  if(any(!is.numeric(arg))) {
      stop(cgMessage("All values of the argument ", argname, " need to be a",
                     if(positive) "greater than zero" else "",
                     "numeric value"))
  }
  if(positive && any(arg <= 0))  {
      stop(cgMessage("All values of the argument ", argname, " need to be of",
                     "positive value."))
  }
  if(integer && any(arg%%1 != 0))  {
      stop(cgMessage("All values of the argument ", argname, " need to be of",
                     "integer value."))
  }
  return(TRUE)
}

validNumericOrCensored <- function(x) {
  ##
  ## PURPOSE: Is every value in the vector coerceable from character to numeric value, or
  ## represents a censored observation via <x or >x ?
  ## Beware that this assumes no missing values and thus does not test for them.
  for(i in seq(along=x)) {
    thevalue <- x[i]
    if(!is.na(as.numeric(thevalue))) return(TRUE)
    else if(regexpr("^\\d+\\.*\\d*$", trimWhiteSpace(thevalue), perl=TRUE))
      return(TRUE)
    else {
      stop(cgMessage("The observation ", x[i],
                     "does not seem to be either coerceable to numeric or",
                     "representative of a censored representation such",
                     "as '<x' or '>x'."))
    }
  }
}


validList <- function(arg, names=NULL, argname=NULL) {
  ##
  ## PURPOSE: Simple check that a list was specified
  ##
  if(is.null(argname)) argname <- deparse(substitute(arg))
  if(!is.list(arg)) {
    stop(cgMessage("The argument ", argname,
                   " needs to be of mode list."))
  }
  if(!is.null(names)) {
    numberofnames <- length(names)
    names.arg <- unique(names(arg))
    matchednames <- 0
    for(i in seq(along=(names.arg))) {
      if(!is.na(pmatch(names.arg[i], names))) {
        matchednames <- matchednames + 1
      }
    }
    if(matchednames < numberofnames) {
      stop(cgMessage("At least one of the names",
                     "needed in the list called", argname,
                     "is not correctly specified or is missing in the function call."))
    }
  }
  return(TRUE)
}

validAtomicVec <- function(arg, lengthminreq=2) {
  ##
  ## PURPOSE: Simple check that a atomic vector of minimum length was specified
  ##
  argname <- deparse(substitute(arg))
  if(!is.atomic(arg) || length(arg) < lengthminreq) {
    stop(cgMessage("The argument ", argname,
                   " needs to be an atomic vector",
                   if(!is.null(lengthminreq)) "of length at least ", lengthminreq))
  }
  else return(TRUE) 
}

validArgMatch <- function(arg, choices, argname=NULL) {
  ##
  ## PURPOSE: Check choices for an argument,
  ## allowing for partial matching
  ##
  argname <- if(!is.null(argname)) argname else deparse(substitute(arg))
  arg <- try(match.arg(arg, choices))
  
  if(class(arg)=="try-error") {
    choicestr <- paste(paste("\"", choices, "\"", sep=""),
                       collapse=" or ")
    choicestr <- paste(choicestr, ".", sep="")
    stop(cgMessage("The", argname, "argument",
                   "needs to be one of", choicestr))
  }
  else return(arg)
}

validDotsArg <- function(dots, arg, choices, default=NULL) {
  ##
  if(!is.na(pmatch(arg, names(dots)))) {
    argname <- arg
    arg <- eval(parse(text=paste("dots$",arg,sep="")))
    arg <- validArgMatch(arg, c(choices, default), argname=argname)
  }
  else {
    arg <- default
  }
  return(arg)
}

validDotsArgs <- function(dots, names) {
  ##
  if(length(dots)==0) return(TRUE)
  dots.names <- names(dots)
  ## If dots is only of length 1 and not named
  if(is.null(dots.names)) {
    stop(cgMessage("There seems to be an unassigned value",
                   "in the dots arguments portion of the function",
                   "call."))
  }
  ## Otherwise
  sapply(dots.names,
         function(x, names) {
           if(is.na(pmatch(x, names))) {
             stop(cgMessage("The", x, "argument",
                            "is not used for this function,",
                            "and can not be in the function call."))
           }
           invisible(TRUE)
         },
         names=names)
  return(TRUE)
}


getDotsArgName <- function(dots, name) {
  dots.names <- names(dots)
  ## Check this for a better way
  ## indx <- try(which(pmatch(dots.names, name)==1), silent=TRUE)
  indx <- which(pmatch(dots.names, name)==1)
  
  thearg <- NA ## initialization
  ## if(!inherits(indx, "try-error") && any(!is.na(indx))) {
  if(any(!is.na(indx))) {
    thearg <- stripmiss(dots.names[indx])
  }

  return(thearg)
}
                          
parsePartialName <- function(names, fullname, prefix=NULL) {
  indx <- try(pmatch(names, fullname), silent=TRUE)
  
  thearg <- parsedarg <- NA ## initialization
  if(!inherits(indx, "try-error") && any(!is.na(indx))) {
    thearg <- stripmiss(names[as.logical(indx)])
    parsedarg <- parse(text=paste(prefix, thearg, sep=""))    
  }
  return(parsedarg)
}

reportInvalidArg <- function(argname) {
  stop(cgMessage("The", argname, "argument",
                 "does not appear to be specified",
                 "with a value as needed in the function call."))
  invisible()
}

validEqualLength <- function(vec1, vec2) {
  if(length(vec1)!=length(vec2)) {
    vec1name <- deparse(substitute(vec1))
    vec2name <- deparse(substitute(vec2))
    stop(cgMessage(paste("The vectors", vec1name, "and", vec2name,
                         "need to be of the same length."), tryAgain))
  }
  return(TRUE)
}


validArgMatch <- function(arg, choices, argname=NULL) {
  ##
  ## PURPOSE: Check choices for an argument,
  ## allowing for partial matching
  ##
  argname <- if(!is.null(argname)) argname else deparse(substitute(arg))
  arg <- try(match.arg(arg, choices))
  
  if(class(arg)=="try-error") {
    choicestr <- paste(paste("\"", choices, "\"", sep=""),
                       collapse=" or ")
    choicestr <- paste(choicestr, ".", sep="")
    stop(cgMessage("The", argname, "argument",
                   "needs to be one of", choicestr))
  }
  else return(arg)
}

validArgDigits <- function(x) {
  if(is.null(x)) return(NULL)
  if(!(is.numeric(x) && is.element(x, 0:4))) {
    stop(cgMessage("The digits argument needs to be",
                   "set to one of NULL, 0, 1, 2, 3,",
                   "or 4."))
  }
  else return(x)
}

validArgModel <- function(...) {
  dotsnames <- names(list(...))
  if(!is.null(dotsnames) && pmatch(dotsnames, "model", nomatch=FALSE) ) {
    stop(cgMessage("The model argument is not relevant",
                   "since it is designed only for",
                   "lm (ols) or rlm (rr) fitted models.")
         )
  }
  return(TRUE)
}

