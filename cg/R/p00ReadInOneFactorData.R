## $Id: p00ReadInOneFactorData.R 6053 2015-02-22 20:23:45Z bpikouni $
## One-Factor Unpaired Groups Case

## Read-in functionality

setClass("cgOneFactorData",
         representation(dfr="data.frame",
                        dfru="data.frame",
                        fmt.dfru="data.frame",
                        has.censored="logical",
                        settings="list"))

prepareCGOneFactorData <- function(dfr, format="listed",
                                   analysisname="", endptname="",
                                   endptunits="",
                                   logscale=TRUE,
                                   zeroscore=NULL,
                                   addconstant=NULL,
                                   rightcensor=NULL,
                                   leftcensor=NULL,
                                   digits=NULL, 
                                   refgrp=NULL, stamps=FALSE) {
  ##
  ## PURPOSE: Read-in the input data frame and convert it to
  ## "univariate style" if needed for the usual S modeling conventions.
  ##
  ## Input argument handling
  dataformat <- validArgMatch(format, c("listed","groupcolumns"))
  validBoolean(logscale)
  validCensor(rightcensor, direction="right")
  validCensor(leftcensor, direction="left")
  validBoolean(stamps)
  digits <- validArgDigits(digits)

  if(!is.data.frame(dfr)) {
    stop(cgMessage("The input data needs to be of class data.frame"))
  }
  if(!is.null(zeroscore) && !is.null(addconstant)) {
    stop(cgMessage("The zeroscore and addconstant arguments",
                   "cannot be specified at the same time.",
                   seeHelpFile("prepareCGOneFactorData")))
  }
  if(logscale==FALSE && (!is.null(zeroscore) || !is.null(addconstant))) {
    stop(cgMessage("Both the zeroscore and addconstant arguments",
                   "must be NULL when logscale=FALSE.",
                   seeHelpFile("prepareCGOneFactorData")))
  }
  ## End input argument handling
 
  if(dataformat=="groupcolumns") {
    validCGOneFacGroupColDfr(dfr)

    ## Transformation to univariate form, called dfru
    grpnames <- names(dfr)
    grp <- rep(grpnames, each=nrow(dfr))

    ## Ensure that factor levels are in order given by the data column headers
    grpf <- factor(grp, levels=grpnames)
    endpt <- unlist(dfr)

    dfru <- data.frame(grpf=grpf, endpt=I(endpt), row.names=NULL)

    ## Remove any NA or "NA" values
    dfru <- dfru[!(is.na(dfru$endpt) | dfru$endpt=="NA" | dfru$endpt==""), ]

    ## Check of all values
    validNumericOrCensored(dfru$endpt)
    
    ## Recognize censored data if specified
    endpt <- dfru[, 2]
    endpt.hasrightcens <- (regexpr(">", endpt) > 0)
    endpt.hasleftcens <- (regexpr("<", endpt) > 0)

    if(sum(endpt.hasrightcens) > 0) rightcensor <- NULL
    if(sum(endpt.hasleftcens) > 0) leftcensor <- NULL

    dfru <- validCGOneFacListedDfr(dfru, rightcensor, leftcensor)    
    
  }

  else if(dataformat=="listed") {
    dfru <- validCGOneFacListedDfr(dfr, rightcensor, leftcensor)
    names(dfru)[1] <- "grpf"
  }

  digits <- as.integer(ifelse(is.null(digits),
                              getNumDigits(stripmiss(dfru[,2])), digits))
    
  has.censored <- if(ncol(dfru) > 2) { TRUE } else FALSE

  fmt.dfru <- dfru
  if(has.censored) {
    endptf <- dfru$endpt
    endpt1f <- dfru$endpt1
    endpt2f <- dfru$endpt2
    status <- dfru$status
    
    ## handle censored representations
    endptf[status==0] <- paste(">", endptf[status==0], sep="")
    endptf[status==2] <- paste("<", endptf[status==2], sep="")
    endptf[status==3] <- paste(endpt1f[status==3],"-",
             endpt2f[status==3], sep="")
    fmt.dfru$endpt <- endptf
    fmt.dfru$endpt1 <- fmt.dfru$endpt2 <- fmt.dfru$status <- NULL
  }
  else {
    names(fmt.dfru) <- names(dfru) <- c("grpf","endpt")
  }

  if(!is.expression(endptname) && endptname=="") { endptname <- "Endpoint" }
    
  ## Drop any unused factor levels
  dfru$grpf <- factorInSeq(as.character(dfru$grpf))
  
  ## Redefine grpnames to make sure any unused factor levels are dropped
  grpnames <- levels(dfru$grpf)
  
  ## Define reference group
  if(!is.null(refgrp)) {
    if(!is.element(refgrp, grpnames)) {
      stop(cgMessage("The refgrp argument does not appear",
                     "to match one of the group labels.",
                     seeHelpFile("prepareCGOneFactorData")))
    }
    ## else keep specified refgrp value
  }
  else {
    refgrp <- grpnames[1]
  }

  ## Do we need to replace 0's with a non-zero score?
  if(logscale && !is.null(zeroscore)) {

    if(has.censored && !is.numeric(zeroscore)) {
        stop(cgMessage("The zeroscore argument can only",
                       "be numeric when there is censoring",
                       "in the data.",
                     seeHelpFile("prepareCGOneFactorData")))
    }
    else if(has.censored && is.numeric(zeroscore)) {
      endpt <- dfru$endpt
      endpt1 <- dfru$endpt1
      endpt2 <- dfru$endpt2
      endptall <- stripmiss(c(endpt, endpt1, endpt2))

      if(all(endptall > 0)) {
        stop(cgMessage("There do not appear to be any zero",
                       "values in the data, so zeroscore should",
                       "be set to NULL.",
                       seeHelpFile("prepareCGOneFactorData")))
      }

      if(zeroscore >= min(endptall[endptall > 0])) {
        stop(cgMessage("The zeroscore must evaluate to a",
                       "value that is less than all the",
                       "the rest of the endpoint values."),
           seeHelpFile("prepareCGOneFactorData"))
      }
      dfru$endpt[dfru$endpt==0] <- zeroscore
      dfru$endpt1[dfru$endpt1==0] <- zeroscore
      dfru$endpt2[dfru$endpt2==0] <- zeroscore

    }

    else { ## no censoring
    
      endpt <- dfru$endpt
    
      if(all(endpt > 0)) {
        stop(cgMessage("There do not appear to be any zero",
                       "values in the data, so zeroscore should",
                       "be set to NULL.",
                       seeHelpFile("prepareCGOneFactorData")))
      }
      zeroscore <- validZeroScore(zeroscore, endpt)
      
      if(zeroscore >= min(endpt[endpt > 0])) {
        stop(cgMessage("The zeroscore must evaluate to a",
                       "value that is less than all the",
                       "the rest of the endpoint values."),
           seeHelpFile("prepareCGOneFactorData"))
      }
      
      dfru$endpt[dfru$endpt==0] <- zeroscore
    }
  }

  ## Or perhaps add a constant to all values so log is defined?
  if(logscale && !is.null(addconstant)) {

    if(has.censored && !is.numeric(addconstant)) {
      stop(cgMessage("The addconstant argument can only",
                     "be numeric when there is censoring",
                     "in the data.",
                     seeHelpFile("prepareCGOneFactorData")))
    }

    else if(has.censored && is.numeric(addconstant)) {
      endpt <- dfru$endpt
      endpt1 <- dfru$endpt1
      endpt2 <- dfru$endpt2
      endptall <- stripmiss(c(endpt, endpt1, endpt2))

      if(all(endptall > 0)) {
        stop(cgMessage("There do not appear to be any zero",
                       "values in the data, so addconstant should",
                       "be set to NULL.",
                       seeHelpFile("prepareCGOneFactorData")))
      }

      dfru$endpt <- dfru$endpt + addconstant
      dfru$endpt1 <- dfru$endpt1 + addconstant
      dfru$endpt2 <- dfru$endpt2 + addconstant

    }

    else { ## no censoring
      endpt <- dfru$endpt
      if(all(endpt > 0)) {
        stop(cgMessage("There do not appear to be any zero",
                       "values in the data, so addconstant should",
                       "be set to NULL.",
                       seeHelpFile("prepareCGOneFactorData")))
      }
      addconstant <- validAddConstant(addconstant, endpt, dfru$grpf)
      dfru$endpt <- endpt + addconstant
    }
  }

  ## Are there only 4 or less unique values in the endpt, suggesting
  ## that Gaussian and equal variance error assumptions are likely to be
  ## seriously violated?
  if(length(unique(dfru$endpt)) < 5) {
    warning(cgMessage("The endpoint data seems to have only 4 or less",
                      "unique values throughout. If you have small sample",
                      "sizes it is likely",
                      "that Gaussian and equal variance assumptions for validity of",
                      "one-factor modeling and analyses will be seriously violated.",
                      warning=TRUE))
  }

  ## If logscale, are all endpoint values greater than 0?
  if(!has.censored && logscale && any(dfru$endpt <= 0)) {
    stop(cgMessage("There seems to be at least one zero or",
                   "negative value in the endpoint data prepared",
                   "for analysis, so a",
                   "log scale analysis cannot be used, as all",
                   "values need to be greater than zero. ",
                   seeHelpFile("prepareCGOneFactorData")))
  }
  else if(has.censored && logscale &&
          (any(stripmiss(dfru$endpt1) <=0) || any(stripmiss(dfru$endpt2) <=0))) {
    stop(cgMessage("There seems to be at least one zero or",
                   "negative value in the endpoint data prepared",
                   "for analysis, so a",
                   "log scale analysis cannot be used, as all",
                   "values need to be greater than zero. ",
                   seeHelpFile("prepareCGOneFactorData")))
  }

  returnObj <- new("cgOneFactorData",
                   dfr=dfr,
                   dfru=dfru,
                   fmt.dfru=fmt.dfru,
                   has.censored=has.censored,   
                   settings=list(
                     analysisname=analysisname,
                     endptname=endptname,
                     endptunits=endptunits,
                     endptscale=if(logscale) { "log" } else { "original" },
                     zeroscore=zeroscore,
                     addconstant=addconstant,
                     rightcensor=rightcensor,
                     leftcensor=leftcensor,
                     digits=digits,
                     grpnames=grpnames,
                     refgrp=refgrp,
                     stamps=stamps))
  
  returnObj
}


validCGOneFacGroupColDfr <- function(dfr) {
  ##
  ## PURPOSE:  Error checking on dfr in groupcolumn format, the input data frame
  ##
  ## Are there at least 2 groups?
  if(ncol(dfr) < 2) {
    stop(cgMessage("There seems to be only one group in the input data.  At",
                   "least two groups are needed."))
  }
  ## Are there duplicate names?
  else if(length(unique(names(dfr))) < ncol(dfr)) {
    stop(cgMessage("The input data seems to have group names that",
                   "are not all distinct."),
         )
  }

  else return(TRUE)
}


validCGOneFacListedDfr <- function(dfr, rightcensor, leftcensor) {
  ##
  ## PURPOSE:  Error checking on dfr in listed format, the input data frame,
  ## and conversion of censored data into 5 column structure for cg.
  if(ncol(dfr) < 2 || ncol(dfr) > 4) {
    stop(cgMessage("The listed input data format needs to be between 2 and 4",
                   "columns. The first column is the group identifier",
                   "the second is the endpoint, and any others are related to",
                   "censored data formats.",
                   seeHelpFile("prepareCGOneFactorData")))
  }
  
  ## Ensure the first column is character and has at least 2 unique values
  if(length(unique(as.character(dfr[,1]))) < 2) {
    stop(cgMessage("The first column of the listed input data format",
                   "needs to have at least two distinct values since",
                   "it is the group identifier. There seems to be only one.",
                   seeHelpFile("prepareCGOneFactorData")))
  }

  ## Ensure the second column is numeric if there is no censoring
  ##if(is.null(leftcensor) && is.null(rightcensor) && !is.numeric(dfr[, 2])) {
  ##  stop(cgMessage("The second column of the listed input data format",
  ##                 "must have all numeric values when there is no censoring.",
  ##                 seeHelpFile("prepareCGOneFactorData")))
  ##}

  ## Handling of censored data formats
  ##
  if(ncol(dfr)==2) {
    ## presumes a format for second column endpoint that may have "<" or ">" for
    ## censored observations, OR
    ## has rightcensor or leftcensor specified
    dfr <- dfr[!is.na(dfr[, 2]), ]
    endpt <- dfr[, 2]

    if(is.numeric(rightcensor)) endpt[endpt==rightcensor] <- paste(">", rightcensor,
                            sep="")
    if(is.numeric(leftcensor)) endpt[endpt==leftcensor] <- paste("<", leftcensor,
                           sep="")
    endpt.hasright <- (regexpr(">", endpt) > 0)
    endpt.hasleft <- (regexpr("<", endpt) > 0)

    ## If no censoring is evident return dfr as unchanged
    if(sum(endpt.hasright) + sum(endpt.hasleft) == 0) return(dfr)
    
    dfr[, 2] <- as.numeric(gsub(">|<","", endpt))
    dfr[, "endpt2"] <- dfr[, "endpt1"] <- dfr[, 2]
    dfr[, "status"] <- rep(1, nrow(dfr))
    
    dfr$status[endpt.hasright] <- 0
    dfr$status[endpt.hasleft] <- 2
    
  }
  
  else if(ncol(dfr)==3) {
    ## remove any observations with critically missing information
    dfr <- dfr[!(is.na(dfr[,2]) | is.na(dfr[,3])), ]
    
    if(is.null(leftcensor) && is.null(rightcensor)) {
      ## must have 0, 1, and 2 each represented
      statuscheck012 <- table(as.numeric(dfr[[3]])) 
      if(length(statuscheck012)!=3 || sum(as.numeric(names(statuscheck012)))!=3) {
        stop(cgMessage("Each of 0, 1, and 2",
                       "must be represented at least once in the",
                       "status variable",
                       "with this listed input data format,",
                       "when BOTH leftcensor=NULL and rightcensor=NULL.",
                       "Perhaps there is only one of these that",
                       "needs to be specified as TRUE?",
                       seeHelpFile("prepareCGOneFactorData")))      
      }
    }

    else {
      ## or else must be a format where there
      ## is left censoring OR right censoring,
      ## but not both
      if(!is.null(leftcensor) && !is.null(rightcensor)) {
        stop(cgMessage("The leftcensor and rightcensor",
                       "arguments cannot both be specified",
                       "with this listed input data format.",
                       "One and only one must be specified as TRUE",
                       "and the other as NULL.",
                       seeHelpFile("prepareCGOneFactorData")))      
      }
              
      statuscheck01 <- table(as.numeric(dfr[[3]])) 
      if(length(statuscheck01)!=2 || sum(as.numeric(names(statuscheck01)))!=1) {
        stop(cgMessage("The censor status third column",
                       "of the listed input data format must contain",
                       "coercable values of 0 or 1, and each",
                       "must occur at least once.",
                       seeHelpFile("prepareCGOneFactorData")))      
      }

    }
      
    status <- dfr[, 3]
    dfr[, 3] <- NULL
    dfr[, "endpt2"] <- dfr[, "endpt1"] <- dfr[, 2]
    dfr[, "status"] <- status

    if(!is.null(leftcensor) && leftcensor) {
      endpt.hasleft <- (status==0)
      dfr[, "status"] <- status
      dfr$status[dfr$status==0] <- 2
    }
  }

  else if(ncol(dfr)==4) {
      
    ## assume order of columns as: grpf, endpt1, endpt2, status
    ## check that the last 3 columns are consistent
    if(!is.null(rightcensor)) {
      warning(cgMessage("The rightcensor argument is ignored",
                        "and set to NULL.", warning=TRUE))
      rightcensor <- NULL
    }
    if(!is.null(leftcensor)) {
      warning(cgMessage("The leftcensor argument is ignored",
                        "and set to NULL.", warning=TRUE))
      leftcensor <- NULL
    }
    

    ## remove any observations with critically missing information
    dfr <- dfr[!(xor((is.na(dfr[,2]) & is.na(dfr[,3])), is.na(dfr[,4]))), ]
    
    thecheck <- function(endpt1, endpt2, status) {
      ## assume one observation
      check.rightcens <- (!is.na(endpt1) & !is.na(endpt2) & status==0)
      check.leftcens <- (!is.na(endpt1) & !is.na(endpt2) & status==2)
      check.intervalcens <- (!is.na(endpt1) & !is.na(endpt2) &
                             (endpt2 > endpt1) & (status==3))
      check.notcens <- (!is.na(endpt1) & !is.na(endpt2) &
                        (endpt2==endpt1) & (status==1))
      if(!(check.rightcens || check.leftcens || check.intervalcens ||
           check.notcens)) {
        FALSE
      }
      else return(TRUE)
    }
    
    for(i in 1:nrow(dfr)) {
      if(!thecheck(dfr[i, 2], dfr[i, 3], dfr[i, 4])) {
        stop(cgMessage("The censored data values for",
                       "Observation / Row", i, "appears to not",
                       "be set up correctly.",
                       seeHelpFile("prepareCGOneFactorData")))      
      }
    }

    
    dfr <- data.frame(grpf=I(dfr[, 1]), endpt=pmin(dfr[, 2], dfr[, 3], na.rm=TRUE),
                      endpt1=dfr[, 2], endpt2=dfr[, 3], status=dfr[, 4])
    dfr$endpt[dfr$status==3] <- NA  ## interval censored too ambiguous to
                                    ## represent with one value
    
  }

  else {
    stop(cgMessage("There is something unexpectedly wrong with the",
                   "listed data structure that cg cannot specifically identify.",
                   seeHelpFile("prepareCGOneFactorData")))      
  }
  
  return(dfr)
}

validCensor <- function(x, direction) {
  if(is.null(x)) { return (TRUE) }
  else if(!(is.numeric(x) || isTRUE(x))) {
    stop(cgMessage(paste("The", paste(direction, "censor", sep=""),
                         "argument needs to be a number, NULL, or TRUE"),
                   seeHelpFile("prepareCGOneFactorData")))
  }
  return(TRUE)
}


validZeroScore <- function(zeroscore, endpt) {
  if(is.null(zeroscore) || is.numeric(zeroscore)) { return(zeroscore) }
  else if(is.character(zeroscore)) {
    if(zeroscore=="estimate") {
      nzendpt <- unique(endpt[endpt > 0])
      ## use extrapolation feature of stats::spline function
      return(exp(stats::spline(nzendpt, log(nzendpt), method="natural",
                        xmin=0, n=length(nzendpt) + 1)$y[1]))
    }
  }
  else {
    stop(cgMessage("The zeroscore argument is not valid.",
                   "If specifying a character value, it needs ",
                   "to be \"estimate\".",
                   "Otherwise, it needs to be NULL or evaluate",
                   "to a positive numeric value that is smaller than",
                   "the rest of the endpoint values."))
  }
  invisible(NULL)
}

validAddConstant <- function(addconstant, endpt, grpf) {
  if(is.null(addconstant) || is.numeric(addconstant)) { return(addconstant) }
  else if(is.character(addconstant)) {
    if(addconstant=="simple") {
      if(min(endpt) > 0) {
        stop(cgMessage("All endpoint values are greater than zero",
                       "so the addconstant value must be set to NULL."))
      }
      ## From S-PLUS white book page 68
      return((max(endpt)-min(endpt))*0.0001)
    }
    else if(addconstant=="VR") {
      if(min(endpt) > 0) {
        stop(cgMessage("All endpoint values are greater than zero",
                       "so the addconstant value must be set to NULL."),
             seeHelpFile("prepareCGOneFactorData"))
      }
      ltf <- MASS::logtrans(endpt ~ grpf, plotit=FALSE)
      indx <- unique(which.max(ltf$y))
      ## y <- endpt ## needed for following evaluation:
      ## alphalength <- length(eval(formals(MASS:::logtrans.default)$alpha))
      alphalength <- length(seq(0.5, 6, by = 0.25))  ## default of MASS::logtrans
      
      if(indx==1 || indx==alphalength) {
        ## expand alpha argument default in MASS::logtrans
        e.alpha <- c(0.001, 0.01, 0.1,
                     0.2, 0.3, 0.4, seq(0.5, 10, seq=0.25)) - min(endpt)
        e.ltf <- MASS::logtrans(endpt ~ grpf,
                                 alpha=e.alpha,
                                 plotit=FALSE)
        e.indx <- which.max(e.ltf$y)
        if(indx==1 || e.indx==length(e.alpha)) {
          warning(cgMessage("Caution: The estimation of addconstant with the VR",
                            "method is on the boundary. This should not be trusted.",
                            "You may want to run",
                            "the MASS::logtrans function independently to",
                            "verify that a reasonable value is chosen,",
                            "or use addconstant=\"simple\" instead",
                            "or use the zeroscore argument approach.",
                            seeHelpFile("prepareCGOneFactorData")))
          return(e.ltf$x[e.indx])
        }
      }
      else {
        return(ltf$x[indx])
      }
    }
  }
  else {
    stop(cgMessage("The addconstant argument is not valid.",
                   "If specifying a character value, it needs ",
                   "to be \"VR\" or \"simple\".",
                   "Otherwise, it needs to be NULL or evaluate",
                   "to a positive numeric value that is larger",
                   "than the absolute value of the minimum."),
         seeHelpFile("prepareCGOneFactorData"))
  }
  invisible(NULL)
}





















