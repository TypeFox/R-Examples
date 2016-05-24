"check.arg.values" <- function(arg.val, options)
#Function to check that arg.val is one of the allowed values
#and to return the position of the argument in the set of values
#that is stored in options
{ kopt <- pmatch(arg.val, options)
  if (is.na(kopt))
     stop("Value for argument, ",arg.val,", is not an allowed option")
  return(kopt)
}

"alldiffs" <- function(predictions, differences = NULL, p.differences = NULL, 
                       sed = NULL, LSD = NULL, backtransforms = NULL, 
                       response = NULL, response.title = NULL, 
                       term = NULL, classify = NULL, 
                       tdf = NULL)
{ #Check arguments
  if (!("predicted.value" %in% colnames(predictions)) || 
      !("standard.error" %in% colnames(predictions)) || !("est.status" %in% colnames(predictions))) 
       warning("Predictions argument does not include the expected column names (e.g. predicted.value)")
  npred <- nrow(predictions)
  if ((!is.null(differences) && !("matrix" %in% class(differences))) ||
      (!is.null(p.differences) && !("matrix" %in% class(p.differences))) || 
      (!is.null(sed) && !("matrix" %in% class(sed))))
       warning("At least one of differences, p.differences and sed is not of type matrix")
  if (!is.null(differences) && !is.null(p.differences) && !is.null(sed))
  { dimens <- c(nrow(differences), nrow(p.differences), nrow(sed), 
                ncol(differences), ncol(p.differences), ncol(sed))
    if (any(npred != dimens))
      stop("At least one of differences, p.differences or sed is not conformable with predictions")
  }
  if (!is.null(backtransforms))
  { if (!("backtransformed.predictions" %in% colnames(backtransforms))) 
       warning("Backtransforms argument does not include a column named backtransformed.predictions")
    if (npred != nrow(backtransforms))
       stop("Backtransforms do not contain the same number of rows as the predictions")
  }
  #ensure diag of sed is NA
  if (!is.null(sed))
    diag(sed) <- NA
  meanLSD <- NULL
  if (!is.null(LSD))
    attr(predictions, which = "meanLSD") <- LSD$meanLSD
  p <- list(predictions = predictions, differences = differences, 
            p.differences = p.differences, sed = sed, LSD = LSD, 
            backtransforms = backtransforms)
  attr(p, which = "response") <- response
  attr(p, which = "response.title") <- response.title
  attr(p, which = "term") <- term
  attr(p, which = "classify") <- classify
  attr(p, which = "tdf") <- tdf
  class(p) <- "alldiffs"
  return(p)
}

"print.alldiffs" <- function(x, which = "all", ...)
 { options <- c("predictions", "backtransforms", 
                "differences", "p.differences", "sed", "LSD", "all")
   opt <- options[unlist(lapply(which, check.arg.values, options=options))]
   title <- attr(x, which = "response.title")
   if (is.null(title))
      title <- as.character(attr(x, which = "response"))
   term <- attr(x, which = "term")
   if (!is.null(term))
      title <- paste(title, " from ", term, sep="")
   #Print predictions and/or LSDs
   if ("all" %in% opt || "predictions" %in% opt || "LSD" %in% opt)
   { if ("all" %in% opt || "predictions" %in% opt)
     { if (!is.null(title))
         cat("\n\n#### Predictions for ", title, "\n\n")
       print(x$predictions)
     } 
     if (!is.null(x$LSD))
     { sed.range <- abs(x$LSD$minLSD - x$LSD$maxLSD) / x$LSD$meanLSD
       cat("\n\nLSD values \n\n")
       cat("minimum LSD = ",x$LSD$minLSD,  "  mean LSD = ",x$LSD$meanLSD,
           "  maximum LSD = ",x$LSD$maxLSD,
           "\n(sed range / mean sed = ",signif(sed.range, digits=3),")\n\n")
     } else
       print(x$LSD)
   }
   if ("all" %in% opt || "differences" %in% opt)
   { cat("\n\nAll pairwise differences between predicted values \n\n")
     print(x$differences, digits=4)
   }
   if ("all" %in% opt || "p.differences" %in% opt)
   { cat("\n\np values for all pairwise differences between predicted values \n\n")
     print(formatC(x$p.differences, digits=3, format="f"), quote=FALSE)
   }
   if ("all" %in% opt || "sed" %in% opt)
   { cat("\n\nStandard errors of differences between predicted values \n\n")
     print(zapsmall(x$sed, 4))
   }
   if (("all" %in% opt & !is.null(x$backtransforms)) || "backtransforms" %in% opt)
   { if (!is.null(title))
        cat("\n\n#### Backtransforms of predictions for ", title, "\n\n")
     print(x$backtransforms)
   }
   invisible()
}

"asrtests" <- function(asreml.obj, wald.tab = NULL, test.summary = NULL, denDF = "default", ...)
{ if (is.null(test.summary))
  { test.summary <- data.frame(matrix(nrow = 0, ncol=5))
    colnames(test.summary) <- c("terms","DF","denDF","p","action")
  }
  else
   if (!is.data.frame(test.summary) || ncol(test.summary) != 5)
     stop("test.summary in an asrtests object should be a data.frame with 5 columns")
  if (!is.null(wald.tab))
  {  if (!is.data.frame(wald.tab) || ncol(wald.tab) != 4)
      stop("wald.tab should be a 4-column data.frame -- perhaps extract Wald component from list")
  }
  else #form wald.tab
  { wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = FALSE, ...)
    if (!is.data.frame(wald.tab))
           wald.tab <- wald.tab$Wald
  }
  if (class(asreml.obj) != "asreml")
    stop("asreml.obj in an asrtests object should be an asreml object")
  test <- list(asreml.obj = asreml.obj, wald.tab=wald.tab, test.summary = test.summary)
  class(test) <- "asrtests"
  return(test)
}


"print.asrtests" <- function(x, which = "all", ...)
 { options <- c("asremlsummary", "pseudoanova", "testsummary", "all")
   opt <- options[unlist(lapply(which, check.arg.values, options=options))]
   if ("asremlsummary" %in% opt | "all" %in% opt)
      print(summary(x$asreml.obj), ...)
   if ("pseudoanova" %in% opt | "all" %in% opt)
   {  cat("\n\n  Pseudo-anova table for fixed terms \n\n")
      print(x$wald.tab, ...)
   }
   if ("testsummary" %in% opt | "all" %in% opt)
   {  cat("\n\n  Sequence of model terms whose status in the model has been investigated \n\n")
      x$test.summary$p <- round(x$test.summary$p, digits=4)
      print(x$test.summary, digits=4, ...)
   }
   invisible()
}

"recalc.wald.tab.asrtests" <- function(asrtests.obj, recalc.wald = FALSE, 
                                     denDF="default", dDF.na = "none", 
                                     dDF.values = NULL, trace = FALSE, ...)
{ if (is.null(asrtests.obj) | class(asrtests.obj) != "asrtests")
    stop("Must supply an asrtests object")
  asreml.obj <- asrtests.obj$asreml.obj
  #Call wald.asreml if recalc.wald is TRUE
  if (recalc.wald)
  { wald.tab <- asreml::wald.asreml(asreml.obj, denDf=denDF, trace = trace, ...)
    if (!is.data.frame(wald.tab))
      wald.tab <- wald.tab$Wald
  }
  else #extract wald.tab from the asrtests object
    wald.tab <- asrtests.obj$wald.tab
  nofixed <- dim(wald.tab)[1]
  hd <- attr(wald.tab, which = "heading")
  options <- c("none", "residual", "maximum", "supplied")
  opt <- options[check.arg.values(dDF.na, options)]
  if (opt == "supplied" & is.null(dDF.values))
    stop('Need to set dDF.values because have set dDF.na = \"supplied\"')
  if (opt == "supplied")
    if (length(dDF.values) != nofixed)
      stop("Number of  dDF.values must be the same as the number of terms in wald.tab")
  den.df <- NA
  if (!("denDF" %in% colnames(wald.tab))) #no denDF
  { wald.tab <- wald.tab[1:(nofixed-1), c(1,3,4)] #Remove Residual line
    wald.tab$F.inc <- wald.tab$'Wald statistic'/wald.tab$Df
    hd[1] <- "Conservative Wald F tests for fixed effects \n"
    attr(wald.tab, which = "heading") <- hd
    
    #Get denDF
    if (opt == "supplied")
    { wald.tab$denDF <- dDF.values
      warning("Supplied dDF.values used for denDF")
    }
    else
      if (opt == "maximum" | opt == "residual") 
      { wald.tab$denDF <- asreml.obj$nedf
        warning("Residual df used for denDF")
      }
  }
  else #have denom. df
  { den.df.na <- is.na(wald.tab$denDF)
    if (any(den.df.na)) #some denDf are NA and so need to use approximate DF
    { if (opt == "supplied")
      { wald.tab$denDF[den.df.na] <- dDF.values[den.df.na]
        warning("At least some supplied dDF.values used for denDF")
      }
      else
      { if (opt == "maximum") 
        { if (!all(den.df.na))
          { den.df <- max(wald.tab$denDF[-1], na.rm=TRUE) 
            warning("Maximum denDF used for some terms")
          }
          else
          { den.df <- asreml.obj$nedf
            warning("Residual df used for denDF for at least some terms")
          }
        }
        else
        { if (opt == "residual")
          { den.df <- asreml.obj$nedf
          warning("Residual df used for denDF for at least some terms")
          }
        }
        wald.tab$denDF[den.df.na] <- den.df
      }
    }
  }
  #Calc Pr
  wald.tab$Pr <- 1 - pf(wald.tab$F.inc, wald.tab$Df, wald.tab$denDF)
  attr(wald.tab, which = "heading") <- hd
  
  return(wald.tab)
}


"as.terms.object" <- function(terms = NULL, asreml.obj = NULL, ...)
{ if (is.character(terms))
  { terms <- as.formula(paste("~ ",terms, sep=""))
  }
  else
  { if (!is.null(terms))  
      terms <- as.formula(terms)
  }
  if (is.null(terms))
    terms.obj <- NULL
  else
    terms.obj <- terms(terms, 
                       keep.order = T, 
                       data = eval(languageEl(asreml.obj$call, which="data"), 
                                   envir = .GlobalEnv), ...)
  return(terms.obj)
}

"rmTermDescription" <- function(term)
{ #Remove description, if there is one, from term in an asreml termlist
  if (length(grep("!", term, fixed=TRUE))!=0) 
    term <- (strsplit(term, "!", fixed=TRUE) )[[1]][1]
  return(term)
}

"findterm" <- function(term, termlist, rmDescription=TRUE)
#This function finds the position of a term in an asreml termlist 
#It strips off stuff to the right of an !, provided it is not to the right of  
#a term starting with R!
{ if (length(termlist) == 0 | is.null(termlist))
  { k <- 0
  }
  else
  { if (rmDescription)
    { if(substr(term, 1, 2) != "R!")  term <- rmTermDescription(term)
      k <- which(sapply(termlist, 
                        FUN=function(kterm, term)
                       { if (substr(kterm, 1, 2) != "R!")  
                            kterm <- rmTermDescription(kterm)
                         haveterm <- setequal(fac.getinTerm(term), 
                                              fac.getinTerm(kterm))
                       }, 
                       term=term))
    }
    else
    { k <- which(sapply(termlist, 
                        FUN=function(kterm, term)
                       { haveterm <- setequal(fac.getinTerm(term), 
                                              fac.getinTerm(kterm))
                       }, 
                       term=term))
    }
    if (length(k) == 0) k <- 0
  }
  return(k)
}

"fac.getinTerm" <- function(term, rmfunction=FALSE)
#function to return the set of factors/variables in a term separated by ':"
{ if (length(term) != 1)
    stop("Multiple terms supplied where only one allowed")
  vars <- unlist(strsplit(term, ":", fixed=TRUE))
  if (rmfunction)
    vars <- unlist(lapply(vars, rmFunction))
  return(vars)
}

"getTermVars" <- function(term)
  #This function gets the vars in each random term from an asreml termlist 
  #It strips off stuff to the right of an !, provided it is not to the right of  
  #a term starting with R!
{ if(substr(term, 1, 2) != "R!")  term <- rmTermDescription(term)
  vars <- fac.getinTerm(term)
  return(vars)
}

"fac.formTerm" <- function(factors)
#function to form a term from a set of factors/variables in a structure
{ term <- paste(factors,collapse=":")
  return(term)
}

"separateFunction" <- function(var)
  #A function to separate the name of a function and the argument to the function
{ #Remove description, if there is one, from term in an asreml termlist
  if (length(grep("(", var, fixed=TRUE))!=0) 
  { var <- (strsplit(var, "(", fixed=TRUE) )[[1]]
    var[2] <- (strsplit(var[2], ")", fixed=TRUE) )[[1]][1]
  }
  return(var)
}

"rmFunction" <- function(var, asreml.obj)
  #A function that returns the variable without any function
{ var <- separateFunction(var)
  if (length(var)==2)
  { var <- var[2]
    #Check for further arguments and strip, if found
    if (length(grep(",", var, fixed=TRUE))!=0) 
    { var <- (strsplit(var, ",", fixed=TRUE) )[[1]]  
      var <- var[1]
    } 
  }  
    return(var)
}

"getVarCode" <- function(var, asreml.obj)
  #A function that returns a code for a variable
  # 0 for an integer or numeric with no function
  # 5 for any other variable with no function
  # Added to these codes is 1,2,3 or 4 for lin, pol, spl and dev, respectively
{ sys.funcs <- c("lin","pol","spl","dev")
  #Does it have a function
  if (length(grep("(", var, fixed=TRUE))!=0)
  { var.comp <- separateFunction(var)
    k <- match(var.comp[1], sys.funcs, nomatch = 0)
    var <- var.comp[2]
  } else
    k <- 0
  type <- class(eval(languageEl(asreml.obj$call, which="data"))[[match(var, 
                                                                       colnames(eval(languageEl(asreml.obj$call, which="data"))))]])
  if (type == "integer" | type == "numeric")
    code <- k
  else
    code <- k + 5
  return(code)
}

"getVarsCodes" <- function(term, asreml.obj)
  #A function to call getVarCode for each of the vars in a term that have been stored in character vector
{ codes <- unlist(lapply(term[1:length(term)], getVarCode, asreml.obj = asreml.obj))
}

"num.recode" <- function(x, new.values)
#function to form a new variate by changing its unique values to new.values
{ x.values <- sort(unique(x))
  nval <- length(x.values)
  if (nval != length(new.values))
  { stop("Must supply a new value for every unique value in the supplied vector")}
  new.x <- x
  for (i in 1:nval)
  { new.x[x == x.values[i]] <- new.values[i]}
  return(new.x)
}

"power.transform" <- function(var.name, power = 1, offset = 0, scale = 1, 
                              titles = NULL, data)
#Function to perform a power transformation on a variable whose name is given as 
#a character string in var.name. The transformed variable is stored in data
{ k <- match(var.name, names(data))
  if (!is.null(titles) & !is.na(match(var.name, names(titles))))
    title <- titles[[var.name]]
  else
    title <- names(data)[k]
  #Get current variable and its name and title
  tvar.name <- var.name
  ttitle <- title
  y <- data[[k]]
  #Apply scale transformation
  if (scale != 1)
  { y <- scale * y
    if (scale == -1)
    { tvar.name <- paste("neg.",tvar.name,sep="")
      ttitle <- paste("Negative of ", ttitle,sep="")
    } else
    { tvar.name <- paste("scaled.",tvar.name,sep="")
      ttitle <- paste("Scaled ", ttitle,sep="")
    }
  }
  #Apply offset
  if (offset != 0)
  { y <- offset + y
    tvar.name <- paste("offset.",tvar.name,sep="")
    if (scale != 1)
      ttitle <- paste("and ", ttitle,sep="")
    ttitle <- paste("Offset  ", ttitle,sep="")
  }
  #Apply power transformation
  if (power != 1)
  {  if (power == 0)
    { tvar.name <- paste("log.",tvar.name,sep="")
      if (min(y, na.rm = "TRUE") < 1e-04)
        stop("Negative values for log transform - could use offset/scale")
      else
        y <- log(y)
      ttitle <- paste("Logarithm of ", ttitle,sep="")
    } else
    { y <- y^power
      if (power == 0.5)
      { tvar.name <- paste("sqrt.",tvar.name,sep="")
        ttitle <- paste("Square root of ", ttitle,sep="")
      } else
        if (power == -1)
        { tvar.name <- paste("recip.",tvar.name,sep="")
          ttitle <- paste("Reciprocal of ", ttitle,sep="")
        }
        else
        { tvar.name <- paste("power.",tvar.name,sep="")
          ttitle <- paste("Power ",power," of ", ttitle, sep="")
        }
    }
    #Add transformed variable to data
    tk <- match(tvar.name, names(data))
    if (is.na(tk))
      tk <- ncol(data) + 1
    data[[tk]]  <- y
    names(data)[tk] <- tvar.name
    #Add transformed title to titles
    nt <- length(titles)
    titles[[nt+1]] <- ttitle
    names(titles)[nt+1] <- tvar.name
  }
  return(list(data = data, tvar.name = tvar.name, titles = titles))
}

"trend.terms.types" <- function(terms, trend.num = NULL, devn.fac = NULL)
#Examines a set of terms to see if, out of devn.fac and trend.num, whether it 
# some terms with devn.fac and/or sum with trend.num.
{ has.devn.fac <- has.trend.num <- FALSE
  if (length(terms) > 0)
  { for (term in terms)
    { factors <- fac.getinTerm(term)
      if (!is.null(devn.fac) && any(devn.fac %in%  factors))
         has.devn.fac <- TRUE
      else
        if (!is.null(trend.num) && any(grepl(trend.num, factors, fixed = TRUE)))
           has.trend.num <- TRUE
    }
  }
  term.types <- list(has.devn.fac = has.devn.fac, has.trend.num = has.trend.num)
  return(term.types)
}

"setvarianceterms.asreml" <- function(call, terms, ignore.suffices = TRUE, constraints = "P", 
                                      initial.values = NA, ...)
  # call is an unevaluated call to asreml (can create using the call function)
  # terms is a vector of characters with the names of the gammas to be set
  # constraints specifies the constraints to be applied to the terms 
  # initial.values specifies the initial values for the terms
  # - ignore.suffices, constraints and initial.values must be of length 1 or 
  #   the same length as terms
  # - if any of constraints or initial.values is set to NA, then they are 
  #   left unchanged for those terms
{ 
  #test for compatibility of arguments
  nt <- length(terms)
  if (length(ignore.suffices) == 1 & nt != 1)
    ignore.suffices <- rep(ignore.suffices, nt)
  if (length(ignore.suffices) != nt)
    stop("ignore.suffices specification is not consistent with terms")
  if (length(constraints) == 1 & nt != 1)
    constraints <- rep(constraints, nt)
  if (length(constraints) != nt)
    stop("constraints specification is not consistent with terms")
  if (length(initial.values) == 1 & nt != 1)
    initial.values <- rep(initial.values, nt)
  if (length(initial.values) != nt)
    stop("initial.values specification is not consistent with terms")

  #add start.values to call and apply constraints to the gammas specified by terms
  start.call <- call
  languageEl(start.call, which = "start.values") <- TRUE
  gamma.start <- eval(start.call, sys.parent())
  gamma.table <- gamma.start$gammas.table
  gammas <- gamma.table$Gamma
  k <- unlist(lapply(1:nt, 
                     FUN=function(i, terms, termslist, ignore.suffices=TRUE)
                     { k <- findterm(terms[i], termslist, rmDescription=ignore.suffices[i])
                       return(k)
                     }, 
                     terms=terms,
                     termslist=gammas, 
                     ignore.suffices=ignore.suffices))
  if (any(k==0))
    stop(paste("Could not find", paste(terms[k==0], collapse=", ")))
  else
  { if (!all(is.na(constraints)))
    { kk <- k[!is.na(constraints)]
      gamma.table$Constraint[kk] <- constraints[!is.na(constraints)]
    }
    if (!all(is.na(initial.values)))
    { kk <- k[!is.na(initial.values)]
      gamma.table$Value[kk] <- initial.values[!is.na(initial.values)]
    }
  }
  #rerun with the unconstrained parameters, adding parameters in ...
  unconst.call <- call
  languageEl(unconst.call, which = "G.param") <- gamma.table
  languageEl(unconst.call, which = "R.param") <- gamma.table
  
  #deal with args coming via ...
  tempcall <- list(...)
  if (length(tempcall)) 
  { for (z in names(tempcall))
    languageEl(unconst.call, which = z) <- tempcall[[z]]
  }
  #Evaluate the call
  new.reml <- eval(unconst.call, sys.parent())
  new.reml$call <- unconst.call
  invisible(new.reml)
}

"newfit.asreml" <- function(asreml.obj, fixed., random., sparse., rcov., 
                           update = TRUE, keep.order = TRUE, 
                           set.terms = NULL, ignore.suffices = TRUE, 
                           constraints = "P", initial.values = NA, ...)
#a function to refit an asreml model with modified model formula
#using either update.asreml or a direct call to asreml
#- the principal difference is that the latter does not enforce the 
#  use of previous values of the variance parameters as initial values.
#- ... is used to pass arguments to asreml.
# set.terms is a vector of characters with the names of the gammas to be set
# constraints specifies the constraints to be applied to the terms 
# initial.values specifies the initial values for the terms
# - ignore.suffices, constraints and initial.values must be of length 1 or 
#   the same length as terms
# - if any of constraints or initial.values is set to NA, then they are 
#   left unchanged for those terms
{ if (is.null(set.terms) & update) #call update
  { asreml.obj <- update(asreml.obj, fixed. = fixed., random. = random., 
                         sparse. = sparse., rcov. = rcov., 
                         keep.order = keep.order, ...)
  }
  else #call asreml with updated formulae - this code modelled on update.asreml
  { "my.update.formula" <- function(old, new, keep.order = TRUE, ...) 
    #function to update a formula
    { env <- environment(as.formula(old))
      tmp <- update.formula(as.formula(old), as.formula(new))
      out <- formula(terms.formula(tmp, simplify = TRUE, keep.order = keep.order))
      environment(out) <- env
      return(out)
    }
    if (is.null(call <- asreml.obj$call) && 
          is.null(call <- attr(asreml.obj, "call"))) 
      stop("Need an object with call component or attribute")
    #Evaluate formulae in case they are stored in a user-supplied object
    if (!is.null(languageEl(call, which = "fixed")))
      languageEl(call, which = "fixed") <- eval(languageEl(call, which = "fixed"))
    if (!is.null(languageEl(call, which = "random")))
      languageEl(call, which = "random") <- eval(languageEl(call, which = "random"))
    if (!is.null(languageEl(call, which = "sparse")))
      languageEl(call, which = "sparse") <- eval(languageEl(call, which = "sparse"))
    if (!is.null(languageEl(call, which = "rcov")))
      languageEl(call, which = "rcov") <- eval(languageEl(call, which = "rcov"))
    #Now update formulae
    if (!missing(fixed.)) 
      languageEl(call, which = "fixed") <- 
                       my.update.formula(as.formula(languageEl(call, which = "fixed")), 
                                         fixed., keep.order = keep.order)
    if (!missing(random.)) 
      languageEl(call, which = "random") <- 
      { if (!is.null(languageEl(call, which = "random"))) 
          my.update.formula(as.formula(languageEl(call, which = "random")), 
                            random., keep.order = keep.order)
        else 
          random.
      }
    if (!missing(sparse.)) 
      languageEl(call, which = "sparse") <- 
      { if (!is.null(languageEl(call, which = "sparse"))) 
           my.update.formula(as.formula(languageEl(call, which = "sparse")), 
                             sparse., keep.order = keep.order)
        else 
          sparse.
      }
    if (!missing(rcov.)) 
      languageEl(call, which = "rcov") <- 
    { if (!is.null(languageEl(call, which = "rcov"))) 
          my.update.formula(as.formula(languageEl(call, which = "rcov")), 
                            rcov., keep.order = keep.order)
      else 
          rcov.
    }

    if(is.null(set.terms))
    { #If R.param and G.param already set, make them NULL
      if (!is.null(languageEl(call, which = "R.param")))
        languageEl(call, which = "R.param") <- NULL
      if (!is.null(languageEl(call, which = "G.param")))
        languageEl(call, which = "G.param") <- NULL
    }
    else
    { #set variance terms
      #test for compatibility of arguments
      nt <- length(set.terms)
      if (length(ignore.suffices) == 1 & nt != 1)
        ignore.suffices <- rep(ignore.suffices, nt)
      if (length(ignore.suffices) != nt)
        stop("ignore.suffices specification is not consistent with set.terms")
      if (length(constraints) == 1 & nt != 1)
        constraints <- rep(constraints, nt)
      if (length(constraints) != nt)
        stop("constraints specification is not consistent with set.terms")
      if (length(initial.values) == 1 & nt != 1)
        initial.values <- rep(initial.values, nt)
      if (length(initial.values) != nt)
        stop("initial.values specification is not consistent with set.terms")
      
      #add start.values to call and unconstrain the gammas specified by set.terms
      start.call <- call
      languageEl(start.call, which = "start.values") <- TRUE
      gamma.start <- eval(start.call, sys.parent())
      gamma.table <- gamma.start$gammas.table 
      gammas <- gamma.table$Gamma
      k <- unlist(lapply(1:nt, 
                         FUN=function(i, set.terms, termslist, ignore.suffices=TRUE)
                         { k <- findterm(set.terms[i], termslist, rmDescription=ignore.suffices[i])
                           return(k)
                         }, 
                         set.terms=set.terms,
                         termslist=gammas, 
                         ignore.suffices=ignore.suffices))
      if (any(k==0))
      { warning(paste("Could not find", paste(set.terms[k==0], collapse=", ")))
        constraints <- constraints[k != 0]
        initial.values <- initial.values[k != 0]
        k <- k[k != 0]
      }
      if (!all(is.na(constraints)))
      { kk <- k[!is.na(constraints)]
        gamma.table$Constraint[kk] <- constraints[!is.na(constraints)]
      }
      if (!all(is.na(initial.values)))
      { kk <- k[!is.na(initial.values)]
        gamma.table$Value[kk] <- initial.values[!is.na(initial.values)]
      }
      #modify call to set varaince parameters
      languageEl(call, which = "G.param") <- gamma.table
      languageEl(call, which = "R.param") <- gamma.table
    }
    
    #deal with args coming via ...
    tempcall <- list(...)
    if (length(tempcall)) 
    { for (z in names(tempcall))
        languageEl(call, which = z) <- tempcall[[z]]
    }
    #Check whether formulae has been reduced to no terms
    if (!is.null(languageEl(call, which = "random")) && 
        length(attr(terms(as.formula(languageEl(call, which = "random"))), "factors")) == 0) 
                 languageEl(call, which = "random") <- NULL
    if (!is.null(languageEl(call, which = "sparse")) && 
          length(attr(terms(as.formula(languageEl(call, which = "sparse"))), "factors")) == 0) 
      languageEl(call, which = "sparse") <- NULL
    if (!is.null(languageEl(call, which = "rcov")) && 
          length(attr(terms(as.formula(languageEl(call, which = "rcov"))), "factors")) == 0) 
      languageEl(call, which = "rcov") <- NULL
    
    #Evaluate the call
    asreml.new.obj <- eval(call, sys.parent())
    asreml.new.obj$call <- call
    #If not converged, issue warning
    if (!asreml.new.obj$converge)
      warning(asreml.new.obj$last.message)
    invisible(asreml.new.obj)
  }
}

"rmboundary.asrtests" <- function(asrtests.obj, trace = FALSE, update = TRUE, 
                                  set.terms = NULL, ignore.suffices = TRUE, 
                                  constraints = "P", initial.values = NA, ...)
#Removes any boundary or singular terms from the fit stored in asreml.obj, 
#one by one from largest to smallest
{ #check input arguments
  if (class(asrtests.obj) != "asrtests")
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
                         ", must supply an asrtests object")
  asreml.obj <- asrtests.obj$asreml.obj
  reml <- asreml.obj$loglik
  test.summary <- asrtests.obj$test.summary
  no.unremoveable <- 0
  #Loop  until have removed all removable boundary terms
  repeat
  { #Find boundary terms
    vcomp <- summary(asreml.obj)$varcomp
    bound.terms <- vcomp$constraint == "Boundary" | 
                   vcomp$constraint == "Singular"
    #                 | comp$constraint == "Fixed"
    #bound.terms[grep("?", vcomp$constraint, fixed=TRUE)] <- TRUE
    vcomp <- vcomp[bound.terms, ]
    nbound <- nrow(vcomp)
    #No boundary terms, so get out
    if (nbound <= 0) break
    #Reduce bound terms to set of bound terms that can be removed
    k <- 1
    while (k <= nbound)
    { term <- rownames(vcomp)[k]
      ranterms.obj <- as.terms.object(languageEl(asreml.obj$call, which="random"), asreml.obj)
      termno <- findterm(term, labels(ranterms.obj))
      if (termno <= 0) #must be an R term or not a recognisable G term
      { vcomp <- vcomp[-k, ]
        k <- k - 1
        nbound <- nbound - 1
        #Store any unremoveable terms
        if (no.unremoveable == 0)
        { terms.unremoveable <- term
          no.unremoveable <- no.unremoveable + 1
        } else
        { if (!(term %in% terms.unremoveable))
          { terms.unremoveable <- c(terms.unremoveable, term)
            no.unremoveable <- no.unremoveable + 1
          }
        }
      }
      k <- k + 1
    }
    #No removeable boundary terms, so get out
    if (nbound <= 0) break
    #If single term, set it as the terms to remove
    if (nbound == 1)
    { term <- rmTermDescription(rownames(vcomp)[1])
    } else #Choose a term to remove
    { #Classify terms as involving random or not because it involves a covariate
      vcomp.vars <- lapply(rownames(vcomp)[1:nbound], getTermVars)
      vcomp.codes <- lapply(vcomp.vars, getVarsCodes, asreml.obj = asreml.obj)
      vcomp <- within(vcomp, {
                          terms.random <- unlist(lapply(vcomp.codes, function(term){all(term == 5)}))
                          varnos <- unlist(lapply(vcomp.codes, length))
                    })
      #Identify terms of the same type that are next to be removed
      max.no.factor <- TRUE
      #Check for terms that involve the dev function
      this.type <- unlist(lapply(vcomp.codes, function(term){any(term == 4 | term == 9)}))
      if (any(this.type))
      { vcomp <- subset(vcomp, this.type)
      } else #Check for random terms (factors only)
      {   this.type <- unlist(lapply(vcomp.codes, function(term){all(term == 5)}))
          if (any(this.type))
          { vcomp <- subset(vcomp, this.type)
            max.no.factor <- FALSE
          } else #Check for terms with spl
          { this.type <- unlist(lapply(vcomp.codes, function(term){any(term == 3 | term == 8)}))
            if (any(this.type))
              vcomp <- subset(vcomp, this.type)
            else #Check for terms with pol
            { this.type <- unlist(lapply(vcomp.codes, function(term){any(term == 2 | term == 7)}))
              if (any(this.type))
                vcomp <- subset(vcomp, this.type)
              else #Check for terms with covariate or lin
              { this.type <- unlist(lapply(vcomp.codes, function(term){any(term < 2 | term == 6)}))
                if (any(this.type))
                  vcomp <- subset(vcomp, this.type)
              }
                
            }
          }
      }
      #Reduce to subset of terms of the type to be removed and choose a term  
      if (max.no.factor)
      { #Get smallest value term amongst those with the most vars
        vcomp <- subset(vcomp, vcomp$varnos == max(vcomp$varnos))
        vcomp <- with(vcomp, vcomp[order(-component),])
        term <- rmTermDescription(tail(rownames(vcomp), 1))
      } else #Get smallest value term amongst those with the least vars
      {   vcomp <- subset(vcomp, vcomp$varnos == min(vcomp$varnos))
          vcomp <- with(vcomp, vcomp[order(-component),])
          term <- rmTermDescription(tail(rownames(vcomp), 1))
      }
    }
    #Remove chosen term
    test.summary <- rbind(test.summary, 
                          data.frame(terms = term, DF = 1, denDF = NA, p = NA, 
                                     action = "Boundary", stringsAsFactors = FALSE))
    mod.ran <- as.formula(paste("~ . - ", term, sep=""))
    asreml.obj <- newfit.asreml(asreml.obj, random. = mod.ran, trace = trace, 
                                update = update, set.terms = set.terms, 
                                ignore.suffices = ignore.suffices, 
                                constraints = constraints, 
                                initial.values = initial.values, ...)
  }
  #Output warning if there are unremoveable bound terms
  if (no.unremoveable > 0)
    warning("\nIn analysing ",asreml.obj$fixed.formula[[2]],
            ", cannot remove the following boundary/singular term(s): ", 
            terms.unremoveable, "\n\n")
  #Check for variance terms that have been fixed and issue warning
  vcomp <- summary(asreml.obj)$varcomp
  if (length(vcomp$constraint[vcomp$constraint == "Fixed"]) != 0)
  { warning("In analysing ",asreml.obj$fixed.formula[[2]],
            ", estimates of the following parameter(s) are fixed:\n",
            rownames(vcomp)[vcomp$constraint == "Fixed"])
  }
  #check for change in log likelihood and if changed issue warning
  change <- abs(reml -asreml.obj$loglik)/reml*100
  if (change > 1)
    warning(paste("Removing boundary terms has changed the log likelihood by "),
            change,"%")
  results <- asrtests(asreml.obj = asreml.obj, 
                      wald.tab = asrtests.obj$wald.tab, 
                      test.summary = test.summary)
  invisible(results)
}

"addrm.terms.asrtests" <- function(terms = NULL, asrtests.obj, add = FALSE, random = FALSE,  
                                   label = NULL, denDF = "default", trace = FALSE, update = TRUE, 
                                   set.terms = NULL, ignore.suffices = TRUE, 
                                   constraints = "P", initial.values = NA, ...)
#Adds or removes a set of terms from either the fixed or random asreml model
{ #check input arguments
  if (is.null(terms))
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
         ", must supply terms to be removed/added")
  else
    if (!is.character(terms))
      stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
           ", must supply terms as character")
  if (class(asrtests.obj) != "asrtests")
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
         ", must supply an asrtests object")
  #initialize
  asreml.obj <- asrtests.obj$asreml.obj
  wald.tab <- asrtests.obj$wald.tab
  test.summary <- asrtests.obj$test.summary
  if (add)
  { init <- " ~ . + ( "
    action <- "Added"
  }
  else
  { init <- " ~ . - ( "
    action <- "Removed"
  }
  #Choose between fixed and random models and update the model
  if (!random)
  { init <- paste(". ", init, sep="")
    term.form <- as.formula(paste(init, terms, " )", sep=""))
    asreml.obj <- newfit.asreml(asreml.obj, fixed. = term.form, 
                                trace = trace, update = update, 
                                set.terms = set.terms, 
                                ignore.suffices = ignore.suffices, 
                                constraints = constraints, 
                                initial.values = initial.values, ...)
    if (is.null(label))
      label <- "Fixed terms"
  }
  else
  { if (is.null(label))
    label <- "Random terms"
    if (is.null(asreml.obj$G.param))
    { if (add)
      { term.form <- as.formula(paste("~ ", terms, sep=""))
        asreml.obj <- newfit.asreml(asreml.obj, random. = term.form, 
                                    trace = trace, update = update, 
                                    set.terms = set.terms, 
                                    ignore.suffices = ignore.suffices, 
                                    constraints = constraints, 
                                    initial.values = initial.values, ...)
      }
      else
      { action <- "Absent"
      }
    }
    else
    { term.form <- as.formula(paste(init, terms, " )", sep=""))
      asreml.obj <- newfit.asreml(asreml.obj, random. = term.form, 
                                  trace = trace, update = update, 
                                  set.terms = set.terms, 
                                  ignore.suffices = ignore.suffices, 
                                  constraints = constraints, 
                                  initial.values = initial.values, ...)
    }
  }
 
  #Update wald.tab
  wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
  if (!is.data.frame(wald.tab))
           wald.tab <- wald.tab$Wald
  #Update test.summary
  if (!asreml.obj$converge)
    action <- paste(action, " - Unconverged", sep="")
  test.summary <- rbind(test.summary, 
                        data.frame(terms = label, DF = NA, denDF = NA, 
                                   p = NA, action = action, 
                                   stringsAsFactors = FALSE))
  results <- asrtests(asreml.obj = asreml.obj, 
                      wald.tab = wald.tab, 
                      test.summary = test.summary)
  invisible(results)
}

"testranfix.asrtests" <- function(term=NULL, asrtests.obj, alpha = 0.05, 
                                  drop.ran.ns = TRUE, positive.zero = FALSE, 
                                  drop.fix.ns = FALSE, denDF="default", dDF.na = "none", 
                                  dDF.values = NULL, trace = FALSE, update = TRUE, 
                                  set.terms = NULL, ignore.suffices = TRUE, 
                                  constraints = "P", initial.values = NA, ...)
#function to test for a single term, using a REMLRT for a random term or based 
#on Wald statistics for a fixed term. Note that fixed terms are never dropped.
{ #check input arguments
  if (class(asrtests.obj) != "asrtests")
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
         ", must supply an asrtests object")
  asreml.obj <- asrtests.obj$asreml.obj
  wald.tab <- asrtests.obj$wald.tab
  test.summary <- asrtests.obj$test.summary
  #Check for multiple terms
  term.form <- as.formula(paste("~ ",term, sep=""))
  term.obj <- as.terms.object(term, asreml.obj)
  if (length(labels(term.obj)) != 1)
   stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
        ", multiple terms not allowed in testranfix.asrtests")
  #Test whether term is in random model
  ranterms.obj <- as.terms.object(languageEl(asreml.obj$call, which="random"), asreml.obj)
  termno <- findterm(term, labels(ranterms.obj))
  if (termno == 0)
  #See if in fixed model
  { termno <- findterm(term, rownames(asreml.obj$aovTbl))
    #Term is not in either model
    if (termno == 0)
      #Term is not in either model
      test.summary <- rbind(test.summary, 
                            data.frame(terms = term, DF = NA, denDF = NA, 
                                       p = NA, action = "Absent", 
                                       stringsAsFactors = FALSE))
    else
    #Have a fixed term
    { wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
      if (!is.data.frame(wald.tab))
             wald.tab <- wald.tab$Wald
      options <- c("none", "residual", "maximum", "supplied")
      opt <- options[check.arg.values(dDF.na, options)]
      if (opt == "supplied" & is.null(dDF.values))
            stop('Need to set dDF.values because have set dDF.na = \"supplied\"')
      #Compute p-value
      p <- wald.tab[termno, 4]
      ndf <- wald.tab$Df[termno]
      den.df <- NA
      if ("denDF" %in% colnames(wald.tab) & !is.na(wald.tab$denDF[termno]))
         den.df <- wald.tab$denDF[termno]
      else
      { if (opt == "supplied")
          den.df <- dDF.values[termno]
        else
        { if (opt == "maximum") 
          { if ("denDF" %in% colnames(wald.tab) & !all(is.na(wald.tab$denDF)))
               den.df <- max(wald.tab$denDF[-1], na.rm=TRUE) 
            else
               den.df <- asreml.obj$nedf
          }
          else
          { if (opt == "residual")
                den.df <- asreml.obj$nedf
            else
                p <- NA
          }
        }
      }
      #Calc F, if necessary, and p
      if (!is.na(den.df))
      { if ("denDF" %in% colnames(wald.tab))
          test.stat <- wald.tab$F.inc[termno]
        else
          test.stat <- wald.tab$'Wald statistic'[termno]/ndf
        p <- 1 - pf(test.stat, ndf, den.df)
      }

      #Add record for test to test.summary and, if drop.fix.ns is TRUE, remove term
      if (is.na(p))
          test.summary <- rbind(test.summary, 
                                data.frame(terms = rownames(wald.tab)[termno], 
                                           DF = ndf, denDF = NA, p = p, 
                                           action = NA, stringsAsFactors = FALSE))
      else
        if (p <= alpha) 
          if (drop.fix.ns)
            test.summary <- rbind(test.summary, 
                                  data.frame(terms = rownames(wald.tab)[termno], 
                                             DF = ndf, denDF = den.df, 
                                             p = p, action = "Retained",
                                             stringsAsFactors = FALSE))
          else
            test.summary <- rbind(test.summary, 
                                  data.frame(terms = rownames(wald.tab)[termno], 
                                             DF = ndf, denDF = den.df, 
                                             p = p, action = "Significant",
                                             stringsAsFactors = FALSE))
      else
        if (drop.fix.ns)
        { test.summary <- rbind(test.summary, 
                                data.frame(terms = rownames(wald.tab)[termno], 
                                           DF = ndf, denDF = den.df, 
                                           p = p, action = "Dropped",
                                           stringsAsFactors = FALSE))
          term.form <- as.formula(paste(". ~ . - ",term, sep=""))
          asreml.obj <- newfit.asreml(asreml.obj, fixed. = term.form, trace = trace, 
                                      update = update, set.terms = set.terms, 
                                      ignore.suffices = ignore.suffices, 
                                      constraints = constraints, 
                                      initial.values = initial.values, ...)
          if (!asreml.obj$converge)
            action <- paste(action, " - Unconverged", sep="")
          #Update wald.tab
          wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
          if (!is.data.frame(wald.tab))
            wald.tab <- wald.tab$Wald
          #Check for boundary terms
          temp.asrt <- rmboundary.asrtests(asrtests(asreml.obj, wald.tab, test.summary), 
                                           trace = trace, update = update, 
                                           set.terms = set.terms, 
                                           ignore.suffices = ignore.suffices, 
                                           constraints = constraints, 
                                           initial.values = initial.values, ...)
          if (nrow(temp.asrt$test.summary) > nrow(test.summary))
            warning("In analysing ",asreml.obj$fixed.formula[[2]],
                    ", Boundary terms removed")
          asreml.obj <- temp.asrt$asreml.obj
          test.summary <- temp.asrt$test.summary
        }
        else
          test.summary <- rbind(test.summary, 
                                data.frame(terms = rownames(wald.tab)[termno], 
                                           DF = ndf, denDF = den.df, 
                                           p = p, action = "Nonsignificant",
                                           stringsAsFactors = FALSE))
    }
  }
  else
  #Remove random term and test
  { term.form <- as.formula(paste("~ . - ",term, sep=""))
    asreml.new.obj <- newfit.asreml(asreml.obj, random. = term.form, trace = trace, 
                                    update = update, set.terms = set.terms, 
                                    ignore.suffices = ignore.suffices, 
                                    constraints = constraints, 
                                    initial.values = initial.values, ...)
    test <- reml.lrt.asreml(asreml.obj, asreml.new.obj, 
                            positive.zero = positive.zero)
    if (test$DF <= 0)
      p <- NA
    else
      p <- test$p
    if (is.na(p) || p <= alpha)
      if (drop.ran.ns)
        test.summary <- rbind(test.summary, 
                            data.frame(terms = term, DF=test$DF, denDF = NA, 
                                       p = p, action = "Retained", 
                                       stringsAsFactors = FALSE))
      else
        test.summary <- rbind(test.summary, 
                            data.frame(terms = term, DF=test$DF, denDF = NA, 
                                       p = p, action = "Significant", 
                                       stringsAsFactors = FALSE))
    else
    { if (drop.ran.ns)
      { test.summary <- rbind(test.summary, 
                            data.frame(terms = term, DF=test$DF, denDF = NA, 
                                       p = p, action = "Dropped", 
                                       stringsAsFactors = FALSE))
        asreml.obj <- asreml.new.obj
        if (!asreml.obj$converge)
           action <- paste(action, " - Unconverged", sep="")
        #Update wald.tab
        wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
        if (!is.data.frame(wald.tab))
           wald.tab <- wald.tab$Wald
        #Check for boundary terms
        temp.asrt <- rmboundary.asrtests(asrtests(asreml.obj, wald.tab, test.summary), 
                                         trace = trace, update = update, 
                                         set.terms = set.terms, 
                                         ignore.suffices = ignore.suffices, 
                                         constraints = constraints, 
                                         initial.values = initial.values, ...)
        if (nrow(temp.asrt$test.summary) > nrow(test.summary))
          warning("In analysing ",asreml.obj$fixed.formula[[2]],
                  ", Boundary terms removed")
        asreml.obj <- temp.asrt$asreml.obj
        test.summary <- temp.asrt$test.summary
    }
      else
        test.summary <- rbind(test.summary, 
                            data.frame(terms = term, DF=test$DF, denDF = NA, 
                                        p = p, action = "Nonsignificant", 
                                       stringsAsFactors = FALSE))
      
    }
  }
  results <- asrtests(asreml.obj = asreml.obj, 
                      wald.tab = wald.tab, 
                      test.summary = test.summary)
  invisible(results)
}

"testswapran.asrtests" <- function(oldterms = NULL, newterms = NULL, asrtests.obj,
                                   label = "Swap in random model", simpler = FALSE, alpha = 0.05, 
                                   positive.zero = FALSE, denDF="default", trace = FALSE, 
                                   update = TRUE, set.terms = NULL, ignore.suffices = TRUE, 
                                   constraints = "P", initial.values = NA, ...)
  #function to test difference between current random model and one in which oldterms are dropped 
  #and newterms are added, using a REMLRT.
{ #check input arguments
  if (class(asrtests.obj) != "asrtests")
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
         ", must supply an asrtests object")
  asreml.obj <- asrtests.obj$asreml.obj
  wald.tab <- asrtests.obj$wald.tab
  test.summary <- asrtests.obj$test.summary

  #Test whether oldterms are in random model
  oldterms.obj <- as.terms.object(oldterms, asreml.obj)
  ranterms.obj <- as.terms.object(languageEl(asreml.obj$call, which="random"), asreml.obj)
  termno <- findterm(oldterms, labels(ranterms.obj))
  if (any(lapply(labels(oldterms.obj), findterm, termlist=labels(ranterms.obj)) == 0))
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
         ", some random terms in oldterms not in random model")
  
  #Remove random oldterms and add random newterms
  new.form <- as.formula(paste("~ . - ",oldterms," + ",newterms, sep=""))
  asreml.new.obj <- newfit.asreml(asreml.obj, random. = new.form, 
                                  trace = trace, update = update, 
                                  set.terms = set.terms, 
                                  ignore.suffices = ignore.suffices, 
                                  constraints = constraints, 
                                  initial.values = initial.values, ...)
  #Perform the test
  change <- FALSE
  if (simpler)
  { test <- reml.lrt.asreml(asreml.obj, asreml.new.obj, 
                            positive.zero = positive.zero)
    if (test$DF <= 0)
      p <- NA
    else
      p <- test$p
    if (!is.na(p) & p <= alpha)
      test.summary <- rbind(test.summary, 
                            data.frame(terms = label, DF=test$DF, denDF = NA, 
                                       p = p, action = "Unswapped", 
                                       stringsAsFactors = FALSE))
    else
    { test.summary <- rbind(test.summary, 
                            data.frame(terms = label, DF=test$DF, denDF = NA, 
                                       p = p, action = "Swapped", 
                                       stringsAsFactors = FALSE))
      change <- TRUE
    }
  }
  else
  { test <- reml.lrt.asreml(asreml.new.obj, asreml.obj, 
                            positive.zero = positive.zero)
    if (test$DF <= 0)
      p <- NA
    else
      p <- test$p
    if (!is.na(p) & p <= alpha)
    { test.summary <- rbind(test.summary, 
                            data.frame(terms = label, DF=test$DF, denDF = NA, 
                                       p = p, action = "Swapped", 
                                       stringsAsFactors = FALSE))
      change <- TRUE
    }
    else
    { test.summary <- rbind(test.summary, 
                            data.frame(terms = label, DF=test$DF, denDF = NA, 
                                       p = p, action = "Rejected", 
                                       stringsAsFactors = FALSE))
    }
  }
  #Update results
  if (change)
  { asreml.obj <- asreml.new.obj
    #check convergence
    if (!asreml.obj$converge)
    { last.act <- length(test.summary$action)
      test.summary$action[last.act] <- paste(test.summary$action[last.act], " - Unconverged", sep="")
    }
    #Check for boundary terms
    temp.asrt <- rmboundary.asrtests(asrtests(asreml.obj, wald.tab, test.summary), 
                                     trace = trace, update = update, 
                                     set.terms = set.terms, 
                                     ignore.suffices = ignore.suffices, 
                                     constraints = constraints, 
                                     initial.values = initial.values, ...)
    if (nrow(temp.asrt$test.summary) > nrow(test.summary))
      warning("In swapping random terms for ",asreml.obj$fixed.formula[[2]],
              ", Boundary terms removed")
    asreml.obj <- temp.asrt$asreml.obj
    test.summary <- temp.asrt$test.summary
    #Update wald.tab
    wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
    if (!is.data.frame(wald.tab))
      wald.tab <- wald.tab$Wald
  }
  results <- asrtests(asreml.obj = asreml.obj, 
                      wald.tab = wald.tab, 
                      test.summary = test.summary)
  invisible(results)
}


"newrcov.asrtests" <- function(terms = NULL, asrtests.obj, label = "R model", 
                               update = TRUE, trace = FALSE, denDF="default", 
                               set.terms = NULL, ignore.suffices = TRUE, 
                               constraints = "P", initial.values = NA, ...)
  #Fits new rcov formula and tests whether the change is significant
{ #check input arguments
  if (is.null(terms))
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
         ", must supply terms to be tested")
  else
    if (!is.character(terms))
      stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
           ", must supply terms as character")
  if (class(asrtests.obj) != "asrtests")
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
         ", must supply an asrtests object")
  #initialize
  asreml.obj <- asrtests.obj$asreml.obj
  wald.tab <- asrtests.obj$wald.tab
  test.summary <- asrtests.obj$test.summary
  term.form <- as.formula(paste("~ ", terms, sep=""))
  #Update the R model
  asreml.obj <- newfit.asreml(asreml.obj, rcov. = term.form, 
                              trace = trace, update = update, 
                              set.terms = set.terms, 
                              ignore.suffices = ignore.suffices, 
                              constraints = constraints, 
                              initial.values = initial.values, ...)
  test.summary <- rbind(test.summary, 
                        data.frame(terms = label, DF=NA, denDF = NA, 
                                   p = NA, action = "Set", 
                                   stringsAsFactors = FALSE))
  #Update results
  #check convergence
  if (!asreml.obj$converge)
  { last.act <- length(test.summary$action)
    test.summary$action[last.act] <- paste(test.summary$action[last.act], " - Unconverged", sep="")
  }
  #Check for boundary terms
  temp.asrt <- rmboundary.asrtests(asrtests(asreml.obj, wald.tab, test.summary), 
                                   trace = trace, update = update, 
                                   set.terms = set.terms, 
                                   ignore.suffices = ignore.suffices, 
                                   constraints = constraints, 
                                   initial.values = initial.values, ...)
  if (nrow(temp.asrt$test.summary) > nrow(test.summary))
    warning("In analysing ",asreml.obj$fixed.formula[[2]],
            ", boundary terms removed")
  asreml.obj <- temp.asrt$asreml.obj
  test.summary <- temp.asrt$test.summary
  #Update wald.tab
  wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
  if (!is.data.frame(wald.tab))
    wald.tab <- wald.tab$Wald
  
    
  results <- asrtests(asreml.obj = asreml.obj, 
                      wald.tab = wald.tab, 
                      test.summary = test.summary)
  invisible(results)
}

"testrcov.asrtests" <- function(terms = NULL, asrtests.obj, label = "R model", 
                                simpler = FALSE, alpha = 0.05, 
                                positive.zero = FALSE, denDF="default", 
                                update = TRUE, trace = FALSE, 
                                set.terms = NULL, ignore.suffices = TRUE, 
                                constraints = "P", initial.values = NA, ...)
#Fits new rcov formula and tests whether the change is significant
{ #check input arguments
  if (is.null(terms))
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
         ", must supply terms to be tested")
  else
    if (!is.character(terms))
      stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
           ", must supply terms as character")
  if (class(asrtests.obj) != "asrtests")
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
         ", must supply an asrtests object")
  #initialize
  asreml.obj <- asrtests.obj$asreml.obj
  wald.tab <- asrtests.obj$wald.tab
  test.summary <- asrtests.obj$test.summary
  term.form <- as.formula(paste("~ ", terms, sep=""))
  #Update the R model
  asreml.new.obj <- newfit.asreml(asreml.obj, rcov. = term.form, 
                                  trace = trace, update = update, 
                                  set.terms = set.terms, 
                                  ignore.suffices = ignore.suffices, 
                                  constraints = constraints, 
                                  initial.values = initial.values, ...)
  change <- FALSE
  if (simpler)
  { test <- reml.lrt.asreml(asreml.obj, asreml.new.obj, 
                            positive.zero = positive.zero)
    if (test$DF <= 0)
      p <- NA
    else
      p <- test$p
    if (!is.na(p) & p <= alpha)
      test.summary <- rbind(test.summary, 
                          data.frame(terms = label, DF=test$DF, denDF = NA, 
                                     p = p, action = "Unswapped", 
                                     stringsAsFactors = FALSE))
    else
    { test.summary <- rbind(test.summary, 
                          data.frame(terms = label, DF=test$DF, denDF = NA, 
                                     p = p, action = "Swapped", 
                                     stringsAsFactors = FALSE))
      change <- TRUE
    }
  }
  else
  { test <- reml.lrt.asreml(asreml.new.obj, asreml.obj, 
                            positive.zero = positive.zero)
    if (test$DF <= 0)
      p <- NA
    else
      p <- test$p
    if (!is.na(p) & p <= alpha)
    { test.summary <- rbind(test.summary, 
                          data.frame(terms = label, DF=test$DF, denDF = NA, 
                                     p = p, action = "Swapped", 
                                     stringsAsFactors = FALSE))
      change <- TRUE
    }
    else
      test.summary <- rbind(test.summary, 
                          data.frame(terms = label, DF=test$DF, denDF = NA, 
                                     p = p, action = "Rejected", 
                                     stringsAsFactors = FALSE))
  }
  #Update results
  if (change)
  { asreml.obj <- asreml.new.obj
    #check convergence
    if (!asreml.obj$converge)
    { last.act <- length(test.summary$action)
      test.summary$action[last.act] <- paste(test.summary$action[last.act], " - Unconverged", sep="")
    }
    #Check for boundary terms
    temp.asrt <- rmboundary.asrtests(asrtests(asreml.obj, wald.tab, test.summary), 
                                     trace = trace, update = update, 
                                     set.terms = set.terms, 
                                     ignore.suffices = ignore.suffices, 
                                     constraints = constraints, 
                                     initial.values = initial.values, ...)
    if (nrow(temp.asrt$test.summary) > nrow(test.summary))
      warning("In analysing ",asreml.obj$fixed.formula[[2]],
              ", boundary terms removed")
    asreml.obj <- temp.asrt$asreml.obj
    test.summary <- temp.asrt$test.summary
    #Update wald.tab
    wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
    if (!is.data.frame(wald.tab))
        wald.tab <- wald.tab$Wald
  }
  results <- asrtests(asreml.obj = asreml.obj, 
                      wald.tab = wald.tab, 
                      test.summary = test.summary)
  invisible(results)
}


"permute.square" <- function(x, permutation)
{ #function to permute the rows and coluns of a square matrix
  if (!is.matrix(x) || nrow(x) != ncol(x))
    stop("x must be a square matrix")
  permuted <- x[permutation, permutation]
  if (!is.null(rownames(x)))
    rownames(x) <- (rownames(x))[permutation]
  if (!is.null(colnames(x)))
    colnames(x) <- (colnames(x))[permutation]
  return(permuted)
}

"permute.to.zero.lowertri" <- function(x)
{ #function to permute a square matrix until all the lower triangular elements are zero
  if (!is.matrix(x) || nrow(x) != ncol(x))
    stop("x must be a square matrix")
  m <- dim(x)[1]
  noperm <- m*(m-1)/2
  if (length(which(x == 0)) < noperm)
    stop("not enough zero elements to permute to all-zero, lower triangular form")
  x.lt <- matrix(0, nrow=m, ncol=m)
  x.lt[lower.tri(x.lt)] <- x [lower.tri(x)]
  i = 0
  #perform permutations until zero lower triangle 
  #or have done twice the no. of possible permutations
  while (any(x.lt == 1) && i <= noperm*2)
  { #find the non-zero element in rightmost column (c) of the lowest row (r)
    nonzero <- which(x.lt == 1, arr.ind=TRUE)
    nonzero <- nonzero[nonzero[, 1] == max(nonzero[, 1]), ]
    if (!is.null(nrow(nonzero)))
        nonzero <- nonzero[nonzero[, 2] == max(nonzero[, 2]), ]
    #perform the permutation that in interchanges rows r and c and columns r and c
    permtn <- 1:m
    permtn[nonzero[1]] <- nonzero[2]
    permtn[nonzero[2]] <- nonzero[1]
    x <- permute.square(x, permtn)
    #get new lower triangle
    x.lt <- matrix(0, nrow=m, ncol=m)
    x.lt[lower.tri(x.lt)] <- x [lower.tri(x)]
  } 
  if (any(x.lt == 1))
    stop("Unable to form matrix with all zero lower triangular elements")
  return(x)
}

"choose.model.asrtests" <- function(terms.marginality=NULL, asrtests.obj, alpha = 0.05, 
                                    drop.ran.ns=TRUE, positive.zero = FALSE, 
                                    drop.fix.ns=FALSE, denDF = "default",  dDF.na = "none", 
                                    dDF.values = NULL, trace = FALSE, update = TRUE, 
                                    set.terms = NULL, ignore.suffices = TRUE, 
                                    constraints = "P", initial.values = NA, ...)
#function to determine the set of significant terms taking into account marginality relations
#terms.marginality should be a square matrix of ones and zeroes with row and column names 
#   being the names of the terms. The diagonal elements should be one, indicating 
#   that a term is marginal to itself. Elements should be one if the row term is 
#   marginal to the column term. All other elements should be zero. 
{ #check input arguments
  if (class(asrtests.obj) != "asrtests")
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
         ", must supply an asrtests object")
  #check matrix  
  if (!is.matrix(terms.marginality) || 
           nrow(terms.marginality) != ncol(terms.marginality))
    stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
         ", terms.marginality must be a square matrix")
  else
   if (is.null(rownames(terms.marginality)) || is.null(colnames(terms.marginality)))
     stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
          ", terms.marginality must have row and column names that are the terms to be tested")
   else
     if (det(terms.marginality) == 0)
       warning("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
               ", Suspect marginalities of terms not properly specified - check")
  #make sure have a terms.marginality matrix with lower triangle all zero
  terms.marginality <- permute.to.zero.lowertri(terms.marginality)
  #perform tests
  sig.terms <- vector("list", length = 0)
  noterms <- dim(terms.marginality)[1]
  current.asrt <- asrtests.obj
  j <- noterms
  #traverse the columns of terms.marginality
  while (j > 0)
  { #get p-value for term for column j and, if random, drop if ns and drop.ran.ns=TRUE
    term <- (rownames(terms.marginality))[j]
    current.asrt <- testranfix.asrtests(term, asrtests.obj = current.asrt, 
                                        alpha=alpha, drop.ran.ns = drop.ran.ns, 
                                        positive.zero = positive.zero, 
                                        drop.fix.ns = drop.fix.ns, 
                                        denDF = denDF, dDF.na = dDF.na, 
                                        dDF.values = dDF.values, trace = trace, 
                                        update = update, set.terms = set.terms, 
                                        ignore.suffices = ignore.suffices, 
                                        constraints = constraints, 
                                        initial.values = initial.values, ...)
    test.summary <- current.asrt$test.summary
    p <- (test.summary[tail(findterm(term, as.character(test.summary$terms)),1), ])$p
    #if significant, add to sig term list and work out which can be tested next
    if (!is.na(p)) 
    { if (p <= alpha)
      { sig.terms <- c(sig.terms, term)
        nonnest <- which(terms.marginality[1:(j-1), j] == 0)
        noterms <- length(nonnest)
        if (noterms == 0)
          j = 0
        else
        { if (noterms == 1)
          { nonnest.name <- rownames(terms.marginality)[nonnest]
            terms.marginality <- terms.marginality[nonnest, nonnest]
            dim(terms.marginality) <- c(1,1)
            rownames(terms.marginality) <- nonnest.name
          }
          else
          { terms.marginality <- terms.marginality[nonnest, nonnest]
          }
          j <- noterms
        }
      }
      #if not significant, proceed to previous column
      else
        j <- j - 1
    }
    else
      j <- j - 1
  }
  invisible(list(asrtests.obj = current.asrt, sig.terms = sig.terms))  
}

"sig.devn.reparam.asrtests" <- function(terms = NULL, asrtests.obj, 
                                        trend.num = NULL, devn.fac = NULL, 
                                        denDF = "default", trace = FALSE, update = TRUE, 
                                        set.terms = NULL, ignore.suffices = TRUE, 
                                        constraints = "P", initial.values = NA, ...)
#reparamterizes a deviations term to a fixed term
#It assumes that the deviations term are deviations from trend in trend.num and 
#  that there is a random deviations term that involves devn.fac
#Also assumes that the same term, but with the devn.fac substitued by trend.num
#   and spl(trend.num), are in the fixed and random models respectively.
#Retains the trend.num term if there are any significant terms involving it.
{ #find out if any terms have either a devn.fac or a trend.num term 
  #  - (marginality implies cannot be both)
  term.types <- trend.terms.types(terms = terms, devn.fac = devn.fac, trend.num = trend.num)
  #Check have terms involving devn.fac or trend.num
  if (term.types$has.devn.fac || term.types$has.trend.num)
    #Now reparamaterize deviations terms, the method depending on whether or not have a mixture
  { for (term in terms)
    { factors <- fac.getinTerm(term)
      #Check this term involves a devn fac
      if (!is.null(devn.fac) && any(devn.fac %in%  factors))
      { #Have mixture so convert random deviations to lin(trend.num) + fixed devn.fac deviations
        if (term.types$has.devn.fac && term.types$has.trend.num)
        { #Add time.num term to be sure - will do nothing if already there
          lin.term <- sub(devn.fac,trend.num,term)
          asrtests.obj <- addrm.terms.asrtests(lin.term, asrtests.obj, add = TRUE, 
                                               random=FALSE, trace = trace, 
                                               update = update, set.terms = set.terms, 
                                               ignore.suffices = ignore.suffices, 
                                               constraints = constraints, 
                                               initial.values = initial.values, ...)
        }
        #remove spl(time.num) and random deviation
        spl.term <- sub(devn.fac,paste("spl(",trend.num,")", sep=""),term)
        ran.term <- paste(spl.term, term, sep = " + " )
        asrtests.obj <- addrm.terms.asrtests(ran.term, asrtests.obj, add=FALSE, 
                                             random=TRUE, trace = trace, 
                                             update = update, set.terms = set.terms, 
                                             ignore.suffices = ignore.suffices, 
                                             constraints = constraints, 
                                             initial.values = initial.values, ...)
        #add devn term to fixed model
        asrtests.obj <- addrm.terms.asrtests(term, asrtests.obj, add=TRUE, 
                                             random=FALSE, denDF = denDF, 
                                             trace = trace, update = update, 
                                             set.terms = set.terms, 
                                             ignore.suffices = ignore.suffices, 
                                             constraints = constraints, 
                                             initial.values = initial.values, ...)
      }
    }
  }
  invisible(asrtests.obj)
}

"predictparallel.asreml" <- function(classify, term = NULL, asreml.obj = NULL, 
                           titles = NULL, x.num = NULL, x.fac = NULL,  
                           x.pred.values = NULL, x.plot.values = NULL, 
                           error.intervals = "Confidence", avsed.tolerance = 0.25, 
                           pairwise = TRUE, tables = "all", levels.length = NA, 
                           transform.power = 1, offset = 0, scale = 1, 
                           inestimable.rm = TRUE, wald.tab = NULL, alpha = 0.05, 
                           dDF.na = "residual", dDF.values = NULL, trace = FALSE, ...)
#a function to get asreml predictions when there a parallel vector and factor are involved
{ if (!is.null(x.pred.values) && !is.null(x.plot.values))
    if (length(x.pred.values) != length(x.plot.values))
       stop("In analysing ",asreml.obj$fixed.formula[[2]],
            ", length of x.pred.values and x.plot.values should be equal")
  if (avsed.tolerance <0 || avsed.tolerance > 1)
    stop("avsed.tolerance should be between 0 and 1")
  int.options <- c("none", "Confidence", "StandardError", "halfLeastSignificant")
  int.opt <- int.options[check.arg.values(error.intervals, int.options)]
  #Get table option and  check if must form pairwise differences
  tab.options <- c("none", "predictions", "backtransforms", 
                   "differences", "p.differences", "sed", "LSD", "all")
  table.opt <- tab.options[unlist(lapply(tables, check.arg.values, options=tab.options))]
  if ("all" %in% table.opt || "differences" %in% table.opt || 
        int.opt == "halfLeastSignificant")
    pairwise <- TRUE
  #Make sure no functions in classify
  vars <- fac.getinTerm(classify, rmfunction=TRUE)
  classify <- fac.formTerm(vars)
  #Get the predicted values when x.num is not involved in classify
  if (is.null(x.num) || !(x.num %in% vars))
  { pred <- predict(asreml.obj, classify=classify, sed=pairwise, trace = trace, ...)$predictions
    if (!is.null(x.fac) && x.fac %in% vars)
    { k <- match(x.fac, names(pred))
      #Set x values for plotting and table labels in x.fac
      if (is.numeric(pred$pvals[[k]]))
      { if (!is.null(x.plot.values))
           pred$pvals[[k]] <- num.recode(pred$pvals[[k]], x.plot.values)
      }
      else
      { if (is.character(pred$pvals[[k]]))
          pred$pvals[[k]] <- factor(pred$pvals[[k]])
        if (!is.null(x.plot.values))
        { pred$pvals[[k]] <- fac.recode(pred$pvals[[k]], x.plot.values)
          pred$pvals[[k]] <- as.numfac(pred$pvals[[k]])
        }
        else
          if (all(is.na((as.numfac(pred$pvals[[k]])))))
            pred$pvals[[k]] <- unclass(pred$pvals[[k]])
          else
            pred$pvals[[k]] <- as.numfac(pred$pvals[[k]])
      }
    }
    else #make sure all variables are factors
    { pred$pvals[vars] <- lapply(1:length(vars), 
                            function(k, vars, data){
                                  if (!is.factor(data[[vars[k]]]))
                                     data[vars[k]] <- factor(data[[vars[k]]])
                                  return(data[[vars[k]]])}, 
                            data=pred$pvals, vars=vars)
    }
  }
  else
  #Get the predicted values when x.num is involved in classify
  { if (is.null(x.pred.values)) 
         warning("In analysing ",asreml.obj$fixed.formula[[2]],
                 ", predictions involve ",x.num," - you may want to specify x.pred.values")
    x.list <- list(x.pred.values)
    names(x.list) <- x.num
    pred <- predict(asreml.obj, classify=classify, levels=x.list, 
                    sed = pairwise, trace = trace, ...)$predictions
    k <- match(x.num, names(pred$pvals))
    #Set x values for plotting and table labels in x.num
    if (!is.null(x.plot.values))
    { pred$pvals[[k]] <- num.recode(pred$pvals[[k]], x.plot.values)
    }
  }
  #Check that x.num and x.fac are conformable
  if ((!is.null(x.num) && x.num %in% names(pred$pvals)) && 
      (!is.null(x.fac) && x.fac %in% names(pred$pvals)))
  { kn <- match(x.num, names(pred$pvals))
    kf <- match(x.fac, names(pred$pvals))  
    if (is.factor(pred$pvals[[kf]]))
       if (length(levels(pred$pvals[[kf]])) != length(unique(pred$pvals[[kn]])))
          stop("In analysing ",asreml.obj$fixed.formula[[2]],
               ", length of ",x.num," and the number of levels of ",x.fac," must be equal")
    else
       if (length(unique(pred$pvals[[kf]])) != length(unique(pred$pvals[[kn]])))
          stop("In analysing ",asreml.obj$fixed.formula[[2]],
               ", length of ",x.num," and the number of levels of ",x.fac," must be equal")
  }
  #Get denominator degrees of freedom
  if (is.null(term))
     term <- classify
#  if (int.opt != "none" & int.opt != "StandardError")
  { if (is.null(wald.tab))
      stop("wald tab needs to be set so that denDF values are available")
    i <- findterm(term, rownames(wald.tab))
    if (i == 0)
    { warning("For ",asreml.obj$fixed.formula[[2]],
              ", ",term," is not a fixed term that has been fitted")
      denom.df <- asreml.obj$nedf
    }
    else
    { #Get options for dDF.na
      options <- c("none", "residual", "maximum", "supplied")
      opt <- options[check.arg.values(dDF.na, options)]
      if (opt == "supplied" & is.null(dDF.values))
        stop('Need to set dDF.values because have set dDF.na = \"supplied\"')
      #Compute denom.df
      denom.df <- NA
      if ("denDF" %in% colnames(wald.tab) & !is.na(wald.tab$denDF[i]))
        denom.df <- wald.tab$denDF[i]
      else
      { if (opt == "supplied")
        denom.df <- dDF.values[i]
        else
        { if (opt == "maximum") 
          { if ("denDF" %in% colnames(wald.tab) & !all(is.na(wald.tab$denDF)))
            denom.df <- max(wald.tab$denDF[-1], na.rm=TRUE) 
            else
              denom.df <- asreml.obj$nedf
          }
          else
          { if (opt == "residual")
            denom.df <- asreml.obj$nedf
          }
        }
      }
    }
  }
  #Set up alldiffs object
  response <- as.character(asreml.obj$fixed.formula[[2]])
  if (!is.null(titles) & !is.na(match(response, names(titles))))
    response.title <- titles[[response]]
  else
    response.title <- response
  diffs <- alldiffs(predictions = pred$pvals, sed = pred$sed, 
                    response = response, response.title =  response.title, 
                    term = term, classify = classify, 
                    tdf = denom.df)
  diffs <- predictiondiffs.asreml(classify=term, alldiffs.obj = diffs, 
                                  x.num = x.num, x.fac = x.fac,  
                                  levels.length = levels.length, 
                                  pairwise = pairwise, 
                                  inestimable.rm = inestimable.rm, 
                                  alpha = alpha)
  #Add lower and upper uncertainty limits
  if (int.opt != "none")
  { revert <- FALSE
    if (int.opt == "StandardError")
      diffs$predictions <- within(diffs$predictions, 
               { lower.StandardError.limit <- 
                        diffs$predictions[["predicted.value"]] - 
                        diffs$predictions[["standard.error"]]
                 upper.StandardError.limit <- 
                        diffs$predictions[["predicted.value"]] + 
                        diffs$predictions[["standard.error"]]
               })
    else
      t.value = qt(1-alpha/2, denom.df)
    if (int.opt == "halfLeastSignificant")
    { sed.range <- abs(diffs$LSD$minLSD - diffs$LSD$maxLSD) /  diffs$LSD$meanLSD
      if (sed.range <= avsed.tolerance)
        diffs$predictions <- within(diffs$predictions, 
               { lower.halfLeastSignificant.limit <- 
                       diffs$predictions[["predicted.value"]] - 0.5 * diffs$LSD$meanLSD
                 upper.halfLeastSignificant.limit <- 
                       diffs$predictions[["predicted.value"]] + 0.5 * diffs$LSD$meanLSD
               })
      else
      { warning("The avsed.tolerance is exceeded - reverting to confidence intervals")
        revert <- TRUE
      }
    } 
    if (int.opt == "Confidence" || revert)
      diffs$predictions <- within(diffs$predictions, 
               { lower.Confidence.limit <- diffs$predictions[["predicted.value"]] - 
                             qt(1-alpha/2, denom.df) * diffs$predictions[["standard.error"]]
                 upper.Confidence.limit <- diffs$predictions[["predicted.value"]] + 
                             qt(1-alpha/2, denom.df) * diffs$predictions[["standard.error"]]
               })
    ks <- match("est.status", names(pred$pvals))
    diffs$predictions <- diffs$predictions[, c(1:(ks-1), (ks+1), (ks+2), ks)]
  }
  #Add backtransforms if there has been a transformation
  backtransforms <- NULL
  if (transform.power != 1 || offset != 0 || scale != 1)
  { backtransforms <- diffs$predictions
    kp <- match("predicted.value", names(backtransforms))
    #Check if LSD used for predictions and so need to compute CIs
    if ((strsplit(names(backtransforms)[kp+2], ".", 
                  fixed=TRUE))[[1]][2] == "halfLeastSignificant")
    { names(backtransforms)[kp+2] <- "lower.Confidence.limit" 
      names(backtransforms)[kp+3] <- "upper.Confidence.limit" 
      backtransforms <- within(backtransforms, 
               { lower.Confidence.limit <- diffs$predictions[["predicted.value"]] - 
                        qt(1-alpha/2, denom.df) * diffs$predictions[["standard.error"]]
                 upper.Confidence.limit <- diffs$predictions[["predicted.value"]] + 
                        qt(1-alpha/2, denom.df) * diffs$predictions[["standard.error"]]
               })
    }
    names(backtransforms)[match("predicted.value", names(backtransforms))] <- 
                                                       "backtransformed.predictions"
    #Backtransform predictions and intervals for power transformation
    if (transform.power == 0)
    { backtransforms$backtransformed.predictions <- 
                                       exp(backtransforms$backtransformed.predictions)
      backtransforms[[kp+2]] <- exp(backtransforms[[kp+2]])
      backtransforms[[kp+3]] <- exp(backtransforms[[kp+3]])
    } else
      if (transform.power != 1)
      { backtransforms$backtransformed.predictions <- 
                        backtransforms$backtransformed.predictions^(1/transform.power)
        backtransforms[[kp+2]] <- backtransforms[[kp+2]]^(1/transform.power)
        backtransforms[[kp+3]] <- backtransforms[[kp+3]]^(1/transform.power)
      } 
#    backtransforms <- backtransforms[, c(1:(kp-1), (kp+5), (kp+2):(kp+4))]
    #Backtransform for offset and scale
    if (offset !=0 || scale != 1)
    { backtransforms$backtransformed.predictions <- 
        backtransforms$backtransformed.predictions/scale - offset
      backtransforms[[kp+2]] <- backtransforms[[kp+2]]/scale - offset
      backtransforms[[kp+3]] <- backtransforms[[kp+3]]/scale - offset
    }
    diffs$backtransforms <- backtransforms
  }
  #Outut tables according to table.opt and save alldiff object for current term
  if (!("none" %in% table.opt))
    print(diffs, which = table.opt)
  return(diffs)
}



"predictionplot.asreml" <- function(classify, y, data, x.num = NULL, x.fac = NULL, 
                                    nonx.fac.order = NULL, colour.scheme = "colour", 
                                    panels = "multiple", graphics.device = NULL,
                                    error.intervals = "Confidence", 
                                    titles = NULL, y.title = NULL, 
                                    filestem = NULL, ...)
#a function to plot asreml predictions and associated statistics
{ #Check options
  scheme.options <- c("black", "colour")
  scheme.opt <- scheme.options[check.arg.values(colour.scheme, scheme.options)]
  panel.options <- c("single", "multiple")
  panel.opt <- panel.options[check.arg.values(panels, panel.options)]
  int.options <- c("none", "Confidence", "StandardError", "halfLeastSignificant")
  int.opt <- int.options[check.arg.values(error.intervals, int.options)]
  intervals <- FALSE
  #work out columns for intervals and their names along with labels for them
  intervals <- TRUE
  klow <- pmatch("lower", names(data))
  kupp <- pmatch("upper", names(data))
  if (klow == 0 || kupp == 0)
    stop("Cannot find a column starting with lower and another with upper to plot intervals")
  ylow <- names(data)[klow]
  yupp <- names(data)[kupp]
  low.parts <- (strsplit(ylow, ".", fixed=TRUE))[[1]]
  upp.parts <- (strsplit(yupp, ".", fixed=TRUE))[[1]]
  if (!setequal(low.parts[-1], upp.parts[-1]))
    stop("Names of columns for lower and upper limits are not consistent")
  if (low.parts[2] == "StandardError")
  { labend <- "+/- standard errors"
    abbrev <- "SE"
  }
  else 
  if (low.parts[2] == "Confidence")
  { labend <- "confidence intervals"    
    abbrev <- "CI"
  }
  else
  if (low.parts[2] == "halfLeastSignificant") 
  { labend <- "+/- half mean LSD"
    abbrev <- "LSI"
  }
  else
    stop("Names for interval limit are not consistent with the types of interval allowed")
  #Do plots
  #Make sure no functions in classify
  vars <- fac.getinTerm(classify, rmfunction = TRUE)
  classify <- fac.formTerm(vars)
  non.x.terms <- setdiff(vars, c(x.fac, x.num))
  n.non.x <- length(non.x.terms)
  if ((length(vars) != n.non.x && n.non.x > 2) || (length(vars) == n.non.x && n.non.x > 3)) 
      stop ("Sorry but plotting for prediction terms with more than 3 variables not implemented")
  #Determine plotting order of non.x.vars
  if (!is.null(nonx.fac.order)) 
  { if (!setequal(nonx.fac.order, non.x.terms))
      stop("The factors nonx.fac.order are not the same as the set in the term being plotted")
    non.x.terms <- nonx.fac.order
    nos.lev <- sapply(1:n.non.x, 
      FUN=function(k, non.x.terms, data){return(length(levels(data[[non.x.terms[k]]])))},
      non.x.terms = non.x.terms, data = data)
  }
  else
  #Default order - sort non.x.vars according to their numbers of levels 
  { nos.lev <- 1
    if(n.non.x >0)
    { nos.lev <- sapply(1:n.non.x, 
                        FUN=function(k, non.x.terms, data)
                               {return(length(levels(data[[non.x.terms[k]]])))},
                        non.x.terms = non.x.terms, data = data)
      non.x.terms <- non.x.terms[order(nos.lev, decreasing = TRUE)]
      nos.lev <- nos.lev[order(nos.lev, decreasing = TRUE)]
    }
  }
  #Append x.vars 
  if (length(vars) != n.non.x)
  { x.terms <- setdiff(vars, non.x.terms)
    vars <- c(non.x.terms, x.terms)
  }
  else
    vars <- non.x.terms
  meanLSD <- attr(data, which = "meanLSD")
  #Plot predicted values
  cbPalette <- rep(c("#CC79A7", "#56B4E9", "#009E73", "#E69F00", "#0072B2", "#D55E00", "#000000"), times=2)
  symb <- rep(c(18,17,15,3,13,8,21,9,3,2,11,1,7,5,10,0), times=10)
  if (!is.null(graphics.device) )
      do.call(graphics.device, list(record = FALSE))
  barplot <- FALSE
  if (length(vars) == length(non.x.terms))
  #Do bar chart
  { barplot <- TRUE
    if (!is.null(titles) & !is.na(match(vars[1], names(titles))))
      x.title <- titles[[vars[1]]]
    else
      x.title <- vars[[1]]
    pred.plot <-  ggplot(data=data, ...) + theme_bw() + 
                    scale_x_discrete(x.title) + scale_y_continuous(y.title)
    if (scheme.opt == "black")
      pred.plot <-  pred.plot + theme_bw() +
                         theme(panel.grid.major = element_line(colour = "grey95", 
                                    size = 0.5), 
                               panel.grid.minor = element_line(colour = "grey98", 
                                    size = 0.5))
    if (n.non.x == 1)
    { pred.plot <-  pred.plot + 
                       aes_string(x = vars[1], y = y) + 
                       scale_fill_manual(values=cbPalette) + 
                       theme(axis.text.x=element_text(angle=90, hjust=1))
      if (scheme.opt == "black")
        pred.plot <-  pred.plot + geom_bar(stat="identity", fill="grey50") 
      else
        pred.plot <-  pred.plot + geom_bar(stat="identity", fill=cbPalette[1])
      if (intervals)
      { if (scheme.opt == "black")
            pred.plot <- pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),   
                                        linetype = "solid", colour = "black") 
        else
            pred.plot <- pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),   
                                        linetype = "solid", colour = cbPalette[5]) 
        pred.plot <- pred.plot + 
                        annotate("text", x=Inf, y=-Inf,  hjust=1, vjust=-0.3, size = 2, 
                                 label = paste("Error bars are ", labend, sep=""))
        if (low.parts[2] == "halfLeastSignificant" && !is.null(meanLSD))
           pred.plot <- pred.plot + 
                         annotate("text", x=-Inf, y=-Inf,  hjust=-0.01, vjust=-0.3, size = 2, 
                                 label = paste("mean LSD = ",signif(meanLSD, digits=3)))
      }
    } else
    if (n.non.x == 2)
    { pred.plot <-  pred.plot + 
                       aes_string(x = vars[1], y = y) + 
                       scale_fill_manual(values=cbPalette) + 
                       facet_grid(as.formula(paste("~ ",vars[2],sep=""))) + 
                       theme(axis.text.x=element_text(angle=90, hjust=1))
      if (scheme.opt == "black")
        pred.plot <-  pred.plot + geom_bar(stat="identity", fill="grey50") 
      else
        pred.plot <-  pred.plot + geom_bar(stat="identity", fill=cbPalette[1])
      if (intervals)
      { non.x.lev <- levels(data[[vars[2]]])
        annot <- data.frame(Inf, -Inf, 
                             factor(non.x.lev[length(non.x.lev)], levels = non.x.lev))
        names(annot) <- c(vars[1], y, vars[2])
        if (scheme.opt == "black")
          pred.plot <- pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),   
                                      linetype = "solid", colour = "black") 
        else
          pred.plot <- pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),   
                                        linetype = "solid", colour = cbPalette[5]) 
        if (nos.lev[2] > 6)
           pred.plot <- pred.plot +
                         theme(strip.text.x = element_text(size = 8))
        pred.plot <- pred.plot + 
                        geom_text(data = annot, label = "Error bars are", 
                                  hjust=1, vjust=-1.3, size = 2) +
                        geom_text(data = annot, label = labend, 
                                  hjust=1, vjust=-0.3, size = 2)
        if (low.parts[2] == "halfLeastSignificant" && !is.null(meanLSD)) 
        { annot <- data.frame(-Inf, -Inf, 
                              factor(non.x.lev[1], levels = non.x.lev))
          names(annot) <- c(vars[1], y, vars[2])
          pred.plot <- pred.plot + 
                        geom_text(data = annot, hjust=-0.01, vjust=-0.3, size = 2, 
                                 label = paste("mean LSD = ",signif(meanLSD, digits=3)))
        }
      }
    } else
    if (n.non.x == 3)
    { pred.plot <-  pred.plot + 
                       aes_string(x = vars[1], y = y) + 
                       scale_fill_manual(values=cbPalette) + 
                       facet_grid(as.formula(paste(vars[3]," ~ ",vars[2],sep=""))) + 
                       theme(axis.text.x=element_text(angle=90, hjust=1))
      if (scheme.opt == "black")
        pred.plot <-  pred.plot + geom_bar(stat="identity", fill="grey50") 
      else
        pred.plot <-  pred.plot + geom_bar(stat="identity", fill=cbPalette[1])
      if (nos.lev[2] > 6 || nos.lev[3] > 6)
         pred.plot <- pred.plot +
                       theme(strip.text.x = element_text(size = 8),
                       strip.text.y = element_text(size = 8, angle = -90))
      if (intervals)
      { non.x.lev1 <- levels(data[[vars[3]]])
        non.x.lev2 <- levels(data[[vars[2]]])
        annot <- data.frame(Inf, -Inf, 
                            factor(non.x.lev1[length(non.x.lev1)], levels = non.x.lev1),
                            factor(non.x.lev2[length(non.x.lev2)], levels = non.x.lev2))
        names(annot) <- c(vars[1], y, vars[c(3,2)])
        if (scheme.opt == "black")
          pred.plot <- pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),   
                                      linetype = "solid", colour = "black") 
        else
          pred.plot <- pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),   
                                        linetype = "solid", colour = cbPalette[5]) 
        pred.plot <- pred.plot + 
                        geom_text(data = annot, label = "Error bars are", 
                                  hjust=1, vjust=-2.1, size = 2) +
                        geom_text(data = annot, label = labend, 
                                  hjust=1, vjust=-1.1, size = 2)
        if (low.parts[2] == "halfLeastSignificant" && !is.null(meanLSD)) 
        { annot <- data.frame(-Inf, -Inf, 
                              factor(non.x.lev1[length(non.x.lev1)], levels = non.x.lev1),
                              factor(non.x.lev2[1], levels = non.x.lev2))
          names(annot) <- c(vars[1], y, vars[c(3,2)])
          pred.plot <- pred.plot + 
                        geom_text(data = annot, hjust=-0.01, vjust=-1.1, size = 2, 
                                 label = paste("mean LSD = ",signif(meanLSD, digits=3)))
        }
      }
    }
  }
  else
  #Do line plot
  { if (intervals)
      int.width <- (max(data[[x.num]]) - min(data[[x.num]]))*0.025
    if (x.num %in% vars)
       x.var <- x.num
    else
       x.var <- x.fac
    if (!is.null(titles) & !is.na(match(x.var,names(titles))))
      x.title <- titles[[x.var]]
    else
      x.title <- x.var
    pred.plot <-  ggplot(data=data, ...) + theme_bw() +
                        scale_x_continuous(x.title) + scale_y_continuous(y.title)
    if (scheme.opt == "black")
      pred.plot <-  pred.plot + theme_bw() +
                         theme(panel.grid.major = element_line(colour = "grey95", 
                                    size = 0.5), 
                               panel.grid.minor = element_line(colour = "grey98", 
                                    size = 0.5))
    #If no non.x.terms ignore single & multiple
    if (n.non.x == 0)
    { pred.plot <-  pred.plot + 
                       aes_string(x = x.var, y = y) + 
                       scale_colour_manual(values=cbPalette) + 
                       scale_shape_manual(values=symb) + 
                       geom_point(shape = symb[7])
      if (scheme.opt == "black")
        pred.plot <-  pred.plot + geom_line(colour = "black") 
      else
        pred.plot <-  pred.plot + geom_line(colour = cbPalette[1])
      if (intervals)
      { if (scheme.opt == "black")
          pred.plot <-  pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),  
                                      linetype = "solid", colour = "black", width=int.width) 
        else
          pred.plot <-  pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),  
                                      linetype = "solid", colour = cbPalette[5], width=int.width) 
        pred.plot <- pred.plot + 
                       annotate("text", x=Inf, y=-Inf,  hjust=1, vjust=-0.3, size = 2, 
                                label =  paste("Error bars are ", labend, sep=""))
       if (low.parts[2] == "halfLeastSignificant" && !is.null(meanLSD)) 
         pred.plot <- pred.plot + 
                       annotate("text", x=-Inf, y=-Inf,  hjust=-0.01, vjust=-0.3, size = 2, 
                               label = paste("mean LSD = ",signif(meanLSD, digits=3)))
      }
    } else
    if (panel.opt == "single")
    #Single with non.x.terms
    {  if (n.non.x == 1)
      { if (scheme.opt == "black")
          pred.plot <-  pred.plot + aes_string(x = x.var, y = y, 
                                      linetype=non.x.terms, shape=non.x.terms)
        else
          pred.plot <-  pred.plot + aes_string(x = x.var, y = y, colour = non.x.terms, 
                                      linetype=non.x.terms, shape=non.x.terms) 
        pred.plot <-  pred.plot + 
                         labs(colour=non.x.terms, linetype=non.x.terms, shape=non.x.terms) +
                        # scale_colour_manual(values=cbPalette) + 
                         scale_shape_manual(values=symb) + 
                         geom_line() + 
                         geom_point()
        if (intervals)
        { pred.plot <- pred.plot + 
                         geom_errorbar(aes_string(ymin=ylow, ymax=yupp),   
                                       linetype = "solid", position = position_dodge(1)) +
                         annotate("text", x=Inf, y=-Inf,  hjust=1, vjust=-0.3, size = 2, 
                                  label =  paste("Error bars are ", labend, sep=""))
          if (low.parts[2] == "halfLeastSignificant" && !is.null(meanLSD)) 
          { pred.plot <- pred.plot + 
                           annotate("text", x=-Inf, y=-Inf,  hjust=-0.01, vjust=-0.3, size = 2, 
                                    label = paste("mean LSD = ",signif(meanLSD, digits=3)))
          }
        }
      }
    }
    else #"multiple"
    { if (n.non.x == 1)
      { pred.plot <- pred.plot + 
                         aes_string(x = x.var, y = y) +
                         scale_colour_manual(values=cbPalette) + 
                         scale_shape_manual(values=symb) + 
                         geom_point(shape = symb[7])  + 
                         facet_wrap(as.formula(paste("~ ",non.x.terms,sep="")))
        if (scheme.opt == "black")
          pred.plot <- pred.plot + geom_line(colour = "black")
        else
          pred.plot <- pred.plot + geom_line(colour = cbPalette[1])
        if (intervals)
        { non.x.lev <- levels(data[[non.x.terms]])
          annot <- data.frame(Inf, -Inf, 
                               factor(non.x.lev[length(non.x.lev)], levels = non.x.lev))
          names(annot) <- c(x.var, y, non.x.terms)
          if (scheme.opt == "black")
            pred.plot <- pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),   
                                        linetype = "solid", colour = "black", width=int.width) 
          else
            pred.plot <- pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),   
                                        linetype = "solid", colour = cbPalette[5], width=int.width) 
          pred.plot <- pred.plot + 
                          geom_text(data = annot, label = "Error bars are", 
                                    hjust=1, vjust=-1.3, size = 2) +
                          geom_text(data = annot, label = labend, 
                                    hjust=1, vjust=-0.3, size = 2)
          if (low.parts[2] == "halfLeastSignificant" && !is.null(meanLSD)) 
          { annot <- data.frame(-Inf, -Inf, 
                                factor(non.x.lev[1], levels = non.x.lev))
            names(annot) <- c(x.var, y, non.x.terms)
            pred.plot <- pred.plot + 
                          geom_text(data = annot, hjust=-0.01, vjust=-0.3, size = 2, 
                                   label = paste("mean LSD = ",signif(meanLSD, digits=3)))
          }
        }
      }
      if (n.non.x == 2)
      { pred.plot <- pred.plot + 
                         aes_string(x = x.var, y = y) +
                         scale_colour_manual(values=cbPalette) + 
                         scale_shape_manual(values=symb) + 
                         geom_point(shape = symb[7]) +
                         facet_grid(as.formula(paste(vars[2]," ~ ",vars[1],sep="")))
        if (scheme.opt == "black")
          pred.plot <- pred.plot + geom_line(colour = "black")
        else
          pred.plot <- pred.plot + geom_line(colour = cbPalette[1])
        if (nos.lev[1] > 6 || nos.lev[2] > 6)
           pred.plot <- pred.plot +
                         theme(strip.text.x = element_text(size = 8),
                         strip.text.y = element_text(size = 8, angle = -90))
        if (intervals)
        { non.x.lev1 <- levels(data[[vars[1]]])
          non.x.lev2 <- levels(data[[vars[2]]])
          annot <- data.frame(Inf, -Inf, 
                               factor(non.x.lev1[length(non.x.lev1)], levels = non.x.lev1),
                               factor(non.x.lev2[length(non.x.lev2)], levels = non.x.lev2))
          names(annot) <- c(x.var, y, non.x.terms)
          if (scheme.opt == "black")
            pred.plot <- pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),   
                                        linetype = "solid", colour = "black", width=int.width) 
          else
            pred.plot <- pred.plot + geom_errorbar(aes_string(ymin=ylow, ymax=yupp),   
                                        linetype = "solid", colour = cbPalette[5], width=int.width) 
          pred.plot <- pred.plot + 
                          geom_text(data = annot, label = "Error bars are", 
                                    hjust=1, vjust=-1.3, size = 2) +
                          geom_text(data = annot, label = labend, 
                                    hjust=1, vjust=-0.3, size = 2)
          if (low.parts[2] == "halfLeastSignificant" && !is.null(meanLSD)) 
          { annot <- data.frame(-Inf, -Inf, 
                                factor(non.x.lev1[1], levels = non.x.lev1),
                                factor(non.x.lev2[length(non.x.lev2)], levels = non.x.lev2))
            names(annot) <- c(x.var, y, non.x.terms)
            pred.plot <- pred.plot + 
                          geom_text(data = annot, hjust=-0.01, vjust=-0.3, size = 2, 
                                   label = paste("mean LSD = ",signif(meanLSD, digits=3)))
          }
        }
      }
    }
  }
  print(pred.plot)
  #Automate saving of plots in files
  if (!is.null(filestem))
  { filename <- paste(filestem,".",gsub(":",".",classify, fixed= TRUE),sep="")
    if (barplot)
      filename <- paste(filename,"Bar", sep=".")
    else
      filename <- paste(filename,"Line", sep=".")
    if (intervals)
      filename <- paste(filename,abbrev,sep=".")
    filename <- paste(filename,".png",sep="")
    savePlot(filename = filename, type = "png")
    cat("\n#### Plot saved in ", filename,"\n")
  }
  invisible()
}

"predictiondiffs.asreml" <- function(classify, alldiffs.obj, 
                                     x.num = NULL, x.fac = NULL,  
                                     levels.length = NA, 
                                     pairwise = TRUE, alpha = 0.05,
                                     inestimable.rm = TRUE)
#a function to get a table of asreml predictions and associated statistics
#  for all pairwise differences
{ #Check alldiffs.obj
  if (!("alldiffs" %in% class(alldiffs.obj)))
    stop("Alldiffs.obj is not of class alldiffs")
  if (is.null(alldiffs.obj$predictions))
    stop("No predictions supplied in alldiffs.obj")
  if (is.null(alldiffs.obj$sed))
    stop(paste("No sed supplied in alldiffs.obj \n",
               "- can obtain using sed=TRUE in predict.asreml"))
  predictions <- alldiffs.obj$predictions
  #Retain only estimable predictions
  which.estim <- (predictions$est.status == "Estimable")
  if (inestimable.rm & sum(which.estim) != nrow(predictions))
  { predictions <- predictions[which.estim, ]
    if (nrow(predictions) == 0)
      warning("There are no estimable predictions")
    #Make sure all factors have only observed levels
    predictions[1:ncol(predictions)] <- 
      lapply(1:ncol(predictions), 
             function(k, data)
             { if (is.factor(data[[k]]))
               data[[k]] <- factor(data[[k]])
               return(data[[k]])
             }, predictions)
    alldiffs.obj$predictions <- predictions
    if (!is.null(alldiffs.obj$sed))
    { if (inestimable.rm)
        alldiffs.obj$sed <- alldiffs.obj$sed[which.estim, which.estim]
      diag(alldiffs.obj$sed) <- NA
    }
    #Reset the other components to NULL
    alldiffs.obj$difference <- NULL
    alldiffs.obj$p.difference <- NULL
    alldiffs.obj$LSD <- NULL
    attr(alldiffs.obj, which = "meanLSD") <- NULL
  }
  
  response <- as.character(attr(alldiffs.obj, which = "response"))
  #Form all pairwise differences, if not present and to be stored
  if (is.null(alldiffs.obj$differences) & pairwise)
  {  #determine factors for row and columns name
    #Make sure no functions in classify
    factors <- fac.getinTerm(classify, rmfunction = TRUE)
    classify <- fac.formTerm(factors)
    nfac <- length(factors)
    #Check all factors in classify are in predictions
    if (length(setdiff (factors, names(predictions))) != 0)
    { if (!is.null(response))
      stop("For ",response,
           ", there are factors in the classify argument that do not have columns in alldiffs.obj$predictions")
      else
        stop("There are factors in the classify argument that do not have columns in alldiffs.obj$predictions")
    }
    #Make sure only one of the numeric and factor that are parallel
    if ((!is.null(x.num) && x.num %in% factors) && (!is.null(x.fac) && x.fac %in% factors))
    { k <- match(x.num, names(predictions))
      predictions <- predictions[, -k]
      nfac <- nfac - 1
    }
    pred.diff <- outer(predictions$predicted.value, predictions$predicted.value, "-")
    #Generate row and column names as the combinations of the levels of factors
    if (nfac > 1)
    { pred.faclist <- vector("list", length=nfac)
      nclassify <- ncol(predictions) - 3
      pred.names <- names(predictions)
      kk <- 0
      for (k in 1:nclassify)
      { if (pred.names[k] %in% factors)
      { kk <- kk + 1
        pred.faclist[[kk]] <- predictions[[k]]
        if (is.numeric(pred.faclist[[kk]]))
          pred.faclist[[kk]] <- factor(pred.faclist[[kk]])
        names(pred.faclist)[kk] <- pred.names[k]
      }
      }
      pred.faclist <- lapply(pred.faclist, FUN=function(ff){if (is.character(levels(ff)) && !is.na(levels.length))
        ff <- factor(ff, labels=substr(levels(ff), start=1, stop=levels.length))
        invisible(ff)})
      pred.lev <- levels(fac.combine(pred.faclist, combine.levels=TRUE))
    }
    else
    { k <- match(factors[[1]], names(predictions))
      pred.fac <- predictions[[k]]
      if (is.numeric(pred.fac))
        pred.fac <- factor(pred.fac)
      pred.lev <- levels(pred.fac)
      if (is.character(pred.lev) && !is.na(levels.length))
        pred.lev <- substr(pred.lev, start=1, stop=levels.length)
    }
    if (nrow(alldiffs.obj$sed) != nrow(pred.diff) | 
          ncol(alldiffs.obj$sed) != ncol(pred.diff))
      stop("The matrix of pairwise differences and sed are not conformable")
    if (ncol(alldiffs.obj$sed) != length(pred.lev) | nrow(alldiffs.obj$sed) != length(pred.lev))
      stop(paste("Dimensions of differences and sed not equal to \n",
                 "the number of observed levels combinations of the factors"))
    dimnames(pred.diff) <- list(pred.lev, pred.lev)
    dimnames(alldiffs.obj$sed) <- list(pred.lev, pred.lev)
  } else
  { if (!pairwise)
    { alldiffs.obj["differences"] <- list(NULL)
      alldiffs.obj["p.differences"] <- list(NULL)
    }
  }
  #Check if tdf available
  denom.df <- attr(alldiffs.obj, which = "tdf")
  if (is.null(denom.df))
    warning(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
                  "- p-values and LSDs not calculated"))
  else
  { #calculate p-values, if not present
    if (is.null(alldiffs.obj$p.differences) & pairwise)
    { p.diff <- abs(pred.diff)/alldiffs.obj$sed
      p.diff <- 2*pt(p.diff, df = denom.df, lower.tail = FALSE)
      alldiffs.obj$differences <- pred.diff
      alldiffs.obj$p.differences <- p.diff
    }
    #calculate LSDs, if not present
    if (is.null(alldiffs.obj$LSD) & pairwise)
    { t.value = qt(1-alpha/2, denom.df)
      minLSD <- t.value * min(alldiffs.obj$sed, na.rm = TRUE)
      maxLSD <- t.value * max(alldiffs.obj$sed, na.rm = TRUE)
      meanLSD <- t.value * mean(alldiffs.obj$sed, na.rm = TRUE)
      alldiffs.obj$LSD <- data.frame(minLSD  = minLSD, 
                                     meanLSD = meanLSD, 
                                     maxLSD = maxLSD)
      attr(alldiffs.obj, which = "meanLSD") <- meanLSD
    }
  }
  return(alldiffs.obj)
}

"pred.present.asreml" <- function(terms, asreml.obj = NULL, 
                                  wald.tab = NULL, dDF.na = "residual", 
                                  dDF.values = NULL, 
                                  x.num = NULL, x.fac = NULL, nonx.fac.order = NULL, 
                                  x.pred.values = NULL, x.plot.values = NULL, 
                                  plots = "predictions", panels = "multiple", 
                                  graphics.device = NULL, 
                                  error.intervals = "Confidence", 
                                  avsed.tolerance = 0.25, titles = NULL, 
                                  colour.scheme = "colour", save.plots = FALSE, 
                                  transform.power = 1, offset = 0, scale = 1, 
                                  pairwise = TRUE, 
                                  tables = "all", levels.length = NA, 
                                  alpha = 0.05, inestimable.rm = TRUE,
                                  trace = FALSE, ...)
#This function forms the predictions for each significant term.
#It then presents either a table or a graph based on the predicted values 
# - the decision is based on whether x.fac or x.num is in the term. 
#Might need to supply additional arguments to predict through "..."
#e.g. present
#Must have levels of x.fac in the order in which they are to be plotted
# - with dates, they should be in the form yyyymmdd
#Probably need to supply x.plot.values if x.fac is to be plotted
{ if (!is.null(x.pred.values) && !is.null(x.plot.values))
    if (length(x.pred.values) != length(x.plot.values))
       stop("In analysing ",asreml.obj$fixed.formula[[2]],
            ", length of x.pred.values and x.plot.values should be equal")
  #Check options
  scheme.options <- c("black", "colour")
  scheme.opt <- scheme.options[check.arg.values(colour.scheme, scheme.options)]
  panel.options <- c("single", "multiple")
  panel.opt <- panel.options[check.arg.values(panels, panel.options)]
  plot.options <- c("none", "predictions", "backtransforms", "both")
  plot.opt <- plot.options[check.arg.values(plots, plot.options)]
  if (avsed.tolerance <0 || avsed.tolerance > 1)
    stop("avsed.tolerance should be between 0 and 1")
  int.options <- c("none", "Confidence", "StandardError", "halfLeastSignificant")
  int.opt <- int.options[check.arg.values(error.intervals, int.options)]
  #Determine if there are any model terms with x.num and any with x.fac
  model.terms <- c(attr(terms(asreml.obj$fixed.formula), which="term.labels"),
             attr(terms(asreml.obj$random.formula), which="term.labels"))
  term.types <- trend.terms.types(terms = model.terms, trend.num = x.num , devn.fac = x.fac)
  diff.list <- vector("list", length = length(terms))
  names(diff.list) <- gsub(":", ".", terms, fixed= TRUE)
  for (term in terms) 
  #Get the predictions
  { #Make sure no functions in classify
    factors <- fac.getinTerm(term, rmfunction = TRUE)
    classify.term <- fac.formTerm(factors)
    #Check if current term involves x.num or x.fac
    if ((!is.null(x.num) && x.num %in% factors) || (!is.null(x.fac) && x.fac %in% factors))
      #Check if the set of model terms is mixed for trend and deviations
    { if (term.types$has.trend.num && term.types$has.devn.fac)
      #Mixed
      { if (!(x.num %in% factors))
          classify.term <- paste(classify.term, x.num, sep=":")
        else
          if (!(x.fac %in% factors))
            classify.term <- paste(classify.term, x.fac, sep=":")
      }
      diffs <- predictparallel.asreml(classify = classify.term, term = term, 
                                     asreml.obj = asreml.obj, titles = titles, 
                                     x.num = x.num, x.fac = x.fac,  
                                     x.pred.values = x.pred.values, 
                                     x.plot.values = x.plot.values, 
                                     error.intervals = error.intervals, 
                                     avsed.tolerance = avsed.tolerance, 
                                     pairwise = pairwise, tables = tables, 
                                     levels.length = levels.length, 
                                     transform.power = transform.power, 
                                     offset = offset, scale = scale, 
                                     inestimable.rm = inestimable.rm, 
                                     wald.tab = wald.tab, dDF.na = dDF.na, 
                                     dDF.values = dDF.values, 
                                     alpha = alpha, trace = trace, ...)
    }
    #Not mixed
    else
      #No need for x.pred.values or to modify term
      diffs <- predictparallel.asreml(classify = classify.term, 
                                      asreml.obj = asreml.obj, titles = titles,   
                                      x.num = x.num, x.fac = x.fac,  
                                      x.plot.values = x.plot.values, 
                                      error.intervals = error.intervals, 
                                      avsed.tolerance = avsed.tolerance, 
                                      pairwise = pairwise, tables = tables, 
                                      levels.length = levels.length, 
                                      transform.power = transform.power, 
                                      offset = offset, scale = scale, 
                                      inestimable.rm = inestimable.rm, 
                                      wald.tab = wald.tab, dDF.na = dDF.na, 
                                      dDF.values = dDF.values, 
                                      alpha = alpha, trace = trace, ...)
    new.name <- gsub(":", ".", term)
    diff.list[[new.name]] <- diffs
    #If plot required, set up y.title
    if (plot.opt != "none")
    { #Set y.title
      y.title <- as.character(attr(diffs, which = "response.title"))
      if (is.null(y.title))
        y.title <- as.character(attr(diffs, which = "response"))
      if (save.plots)
        filestem <- as.character(attr(diffs, which = "response"))
      else
        filestem <- NULL
    }
    #Plot predictions if a plot is required and data are not transformed or if 
    #plot of predictions requested
    transformed <- transform.power != 1 | offset != 0 | scale != 1
    if ((!transformed & plot.opt != "none") | plot.opt == "predictions" 
                                            | plot.opt == "both")
    #Plot predictions
    { predictionplot.asreml(classify = classify.term, y = "predicted.value", 
                            data = diffs$predictions, 
                            x.num = x.num, x.fac = x.fac, 
                            nonx.fac.order = nonx.fac.order, colour.scheme = colour.scheme, 
                            panels = panel.opt, graphics.device = graphics.device,
                            error.intervals = error.intervals,  
                            titles = titles, y.title = y.title, 
                            filestem = filestem, ...)
    }
    if (transformed & (plot.opt == "backtransforms"  | plot.opt == "both"))
    #Plot backtransforms
    { bty.title <- paste("Backtransform of ",y.title,sep="")
      if (save.plots)
        filestem <- paste(filestem,".back",sep="")
      else
        filestem <- NULL
      predictionplot.asreml(classify = classify.term, y = "backtransformed.predictions", 
                            data = diffs$backtransforms, 
                            x.num = x.num, x.fac = x.fac,   
                            nonx.fac.order = nonx.fac.order, colour.scheme = colour.scheme,  
                            panels = panel.opt, graphics.device = graphics.device,
                            error.intervals = error.intervals,
                            titles = titles, y.title = bty.title, 
                            filestem = filestem, ...)
    }
  }
  invisible(diff.list)
}


addrm.terms.asreml <- function(...)
{ .Deprecated(new = "addrm.terms.asrtests", package = "asremlPlus")
  invisible()
}
choose.model.asreml <- function(...)
{ .Deprecated(new = "addrm.terms.asrtests", package = "asremlPlus")
  invisible()
}
recalc.wald.tab.asreml <- function(...)
{ .Deprecated(new = "choose.model.asrtests", package = "asremlPlus")
  invisible()
}
rmboundary.asreml <- function( ...)
{ .Deprecated(new = "rmboundary.asrtests", package = "asremlPlus")
  invisible()
}
sig.devn.reparam.asreml <- function(...)
{ .Deprecated(new = "sig.devn.reparam..asrtests", package = "asremlPlus")
  invisible()
}
testranfix.asreml <- function(...)
{ .Deprecated(new = "testranfix.asrtests", package = "asremlPlus")
  invisible()
}
testrcov.asreml <- function(...)
{ .Deprecated(new = "testrcov.asrtests", package = "asremlPlus")
  invisible()
}
testswapran.asreml <- function(...)
{ .Deprecated(new = "testswapran.asrtests", package = "asremlPlus")
  invisible()
}

#Calls for prediction routines:
# - pred.present.asreml calls predict.parallel.asreml & predictionplot.asreml
# - predict.parallel.asreml calls alldiffs & predictiondiffs.asreml