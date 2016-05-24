#' @title Make a \code{data.table} of Tabulated, Aggregated Values and Weights
#' @description An internal function that aggregates a table
#' and merges in weights.
#' @param data DF/DT; passed to \code{envir} in \code{eval}
#' @param values values to tabulate. Anything \code{evalPopArg} can evaluate.
#' @param print variables to tabulate by and include in \code{prVars} in attributes
#' @param adjust variables to tabulate by and include in \code{adVars} in attributes
#' @param formula a formula such as \code{fot ~ sex} or \code{Surv(fot, lex.Xst) ~ sex}
#' @param Surv.response logical, if \code{TRUE} throws error if response in
#' \code{formula} is not a \code{Surv} object and vice versa
#' @param by.other other variables to tabulate by and include 
#' in \code{boVars} in attributes
#' @param custom.levels a named list of values. When "inflating" the data
#' in the cross-join / cartesian join sense (try e.g. \code{merge(1:5, 1:2)}),
#' one can supply the levels to inflate by using this to ensure inflation is full.
#' E.g. data might only have levels present to do inflation analogous to
#' \code{merge(2:5, 1:2)} although \code{merge(1:5, 1:2)} is intended and 
#' needed. 
#' @param custom.levels.cut.low a character string vector of variable names.
#' These variables mentioned in \code{custom.levels} and existing in data
#' or first modified (in data) using \code{cutLow()} (essentially 
#' \code{cut()} with \code{right = FALSE} and returning the lower bounds
#' as values). Handy for aggregating data e.g. to survival intervals.
#' \strong{NOTE}: the appropriate elements in \code{custom.levels} for these 
#' variables must exceptionally contain an extra value as the roof used in
#' cutting, which will not be used in "inflating" the table using a merge.
#' See Examples.
#' @param weights a named list or long-form data.frame of weights. See Examples.
#' @param internal.weights.values the variable to use to compute internal
#' weights; only used if \code{weights = "internal"}.
#' @param enclos the enclosing environment passed on to \code{eval}. Variables
#' not found in \code{data} or searched for here.
#' @param NA.text a character string to display in a \code{warning}
#' if there are any rows with missing \code{values} or \code{adjust} values.
#' \strong{special:} key phrase \code{\%\%NA_COUNT\%\%} in text is replaced
#' with the count of missing observations.
#' E.g. \code{"Missing \%\%NA_COUNTS\%\% observations due to derpness."}
#' @examples 
#' library(survival)
#' 
#' makeWeightsDT <- popEpi:::makeWeightsDT ## this avoids errors during tests
#' 
#' sire <- copy(popEpi::sire)
#' set.seed(1L)
#' sire$sex <- rbinom(nrow(sire), 1, 0.5)
#' ag <- lexpand(sire, birth = "bi_date", entry = "dg_date", exit = "ex_date",
#'               status = status %in% 1:2, pophaz = popmort, pp = FALSE,
#'               aggre = list(sex, agegr = cut(dg_age, c(0,50,75,Inf)), fot), 
#'               fot = seq(0, 5, 1/12))
#' ps <- quote(list(sex, fot))
#' as <- quote(list(agegr))
#' vs <- list(quote(list(pyrs, at.risk)))
#' ws <- list(agegr = c(0.2,0.4,0.4))
#' 
#' #### custom.levels usage
#' fb <- seq(0, 5-1/12, 1/12) ## exclude 5 as no row has that value
#' ag2 <- ag[fot > 0.5,]
#' # repeats fot intervals < 0.5 as empty rows
#' # may be the safest way to do this
#' dt <- makeWeightsDT(ag2, print = ps, adjust = as, 
#'                     values = vs, weights = ws,
#'                     custom.levels = list(fot = fb))
#' ## aggregate from intervals seq(0, 5, 1/12) to 0:5
#' fb2 <- 0:5 ## (this time we include 5 as the roof)       
#' dt <- makeWeightsDT(ag2, print = ps, adjust = as, 
#'                     values = vs, weights = ws,
#'                     custom.levels = list(fot = fb2),
#'                     custom.levels.cut.low = "fot")              
#'                     
#' 
#' #### use of enclos
#' TF <- environment()
#' gender <- factor(ag$sex)
#' dt <- makeWeightsDT(ag, print = quote(gender), adjust = as, 
#'                     values = vs, weights = ws, enclos = TF)
#' ## or NULL: uses calling frame by default.
#' dt <- makeWeightsDT(ag, print = quote(gender), adjust = as, 
#'                     values = vs, weights = ws,
#'                     enclos = NULL)
#' ## passing parent.fram(1) is the same thing (as below),
#' ## but won't pass in testing these examples somehow (but work in real life)
#' # dt <- makeWeightsDT(ag, print = quote(gender), adjust = as, 
#' #                     values = vs, weights = ws,
#' #                     enclos = NULL)                  
#' 
#' #### formula usage
#' form <- Surv(fot, factor(from0to1))~gender
#' dt <- makeWeightsDT(ag, formula = form, Surv.response = TRUE,
#'                     adjust = as, values = vs, weights = ws,
#'                     enclos = NULL)
#'                     
#' ## or
#' form <- Surv(fot, factor(from0to1))~gender + adjust(agegr)
#' dt <- makeWeightsDT(ag, formula = form, Surv.response = TRUE,
#'                     adjust = NULL, values = vs, weights = ws,
#'                     enclos = NULL)
#'                     
#' ## or   
#' form <- from0to1 ~ fot + gender + adjust(agegr)
#' dt <- makeWeightsDT(ag, formula = form, Surv.response = FALSE,
#'                     adjust = NULL, values = vs, weights = ws,
#'                     enclos = NULL)            
#' 
#' form <- from0to1 ~ fot + adjust(agegr) + adjust(sex)
#' ws2 <- list(agegr = c(0.33, 0.33, 0.33), sex = c(0.5, 0.5))
#' dt <- makeWeightsDT(ag, formula = form, Surv.response = FALSE,
#'                     adjust = NULL, values = vs, weights = ws2,
#'                     enclos = NULL)
#' 
#' ## international standard pops
#' ag <- lexpand(sire, birth = "bi_date", entry = "dg_date", exit = "ex_date",
#'               status = status %in% 1:2, pophaz = popmort, pp = FALSE,
#'               aggre = list(sex, agegr = cut(dg_age, c(seq(0, 85, 5), Inf)), fot), 
#'               fot = seq(0, 5, 1/12))
#'               
#' form <- from0to1 ~ fot + adjust(agegr)
#' dt <- makeWeightsDT(ag, formula = form, Surv.response = FALSE,
#'                     adjust = NULL, values = vs, weights = "world_1966_18of5",
#'                     enclos = NULL)
#'                     
#' form <- from0to1 ~ fot + adjust(agegr, sex)
#' dt <- makeWeightsDT(ag, formula = form, Surv.response = FALSE,
#'                     adjust = NULL, values = vs, 
#'                     weights = list(agegr = "nordic_2000_18of5", sex=c(1,1)),
#'                     enclos = NULL)
makeWeightsDT <- function(data, values = NULL, 
                          print = NULL, adjust = NULL,
                          formula = NULL, Surv.response = TRUE,
                          by.other = NULL, custom.levels = NULL, 
                          custom.levels.cut.low = NULL, weights = NULL, 
                          internal.weights.values = NULL, 
                          enclos = NULL, NA.text = NULL) {
  
  # environmentalism -----------------------------------------------------------
  TF <- environment()
  PF <- parent.frame(1L)
  if (missing(enclos) || is.null(enclos)) {
    enclos <- PF 
  } 
  
  enclos <- eval(enclos, envir = TF)
  
  if (!is.environment(enclos)) {
    stop("Argument 'enclos' is not an environment. (Probably internal error, ",
         "meaning you should complain to the package maintainer if you are not",
         "doing something silly.)")
  }
  
  THIS_CALL <- match.call()
  
  ## dataism -------------------------------------------------------------------
  if (!is.data.frame(data)) stop("data must be a data.frame")
  ## tmpDum for convenience. will be deleted in the end. (if no tabulating vars)
  origData <- data
  tmpDum <- makeTempVarName(origData, pre = "dummy_")
  data <- data.table(rep(1L, nrow(origData)))
  setnames(data, 1, tmpDum)

  
  # formula: vars to print and adjust by ---------------------------------------
  
  adSub <- adjust
  adjust <- evalPopArg(data = origData, arg = adSub, DT = TRUE, enclos = enclos)
  
  if (!is.null(formula)) {
    
    foList <- usePopFormula(formula, adjust = adjust, 
                            data = origData, enclos = enclos, 
                            Surv.response = Surv.response)
    print <- foList$print
    adjust <- foList$adjust
  } else {
    
    prSub <- substitute(print)
    print <- evalPopArg(data = origData, arg = prSub, DT = TRUE, enclos = enclos)
    
    
  }
  
  if (length(weights) && length(adjust)) {
    checkWeights(weights, adjust = adjust)
  }
  
  # variables to print by ----------------------------------------------------
  prVars <- tmpDum
  if (length(print) > 0) {
    prVars <- names(print)
    data[, c(prVars) := TF$print]
    data[, c(tmpDum) := NULL]
  } 
  rm(print)
  
  # standardization ----------------------------------------------------------
  ## note: adjust evaluated above with formula
  adVars <- NULL
  if (length(adjust) > 0) {
    adVars <- names(adjust)
    data[, c(adVars) := TF$adjust]
  } 
  rm(adjust)
  
  if (is.null(weights) && length(adVars)) {
    stop("Variables to adjust by were defined but no weights were supplied.")
  }
  
  # variables to sum -----------------------------------------------------------
  if (!is.list(values)) stop("Argument 'values' must be a list ",
                             "(internal error: complain to the package",
                             "maintainer if you see this)")
  values <- lapply(values, function(x) {
    evalPopArg(data = origData, arg = x, DT = TRUE, enclos = enclos)
  })
  for (dt in setdiff(seq_along(values), 1L)) {
    values[[1L]] <- cbind(values[[1L]], values[[dt]])
  }
  values <- values[[1L]]
  vaVars <- NULL
  if (nrow(values) != nrow(data)) {
    stop("mismatch in numers of rows in data (", nrow(data), 
         ") and 'values' (", nrow(values), "). If you see this message, ",
         "complain to the package maintainer.")
  }
  
  if (length(values) > 0) {
    vaVars <- names(values)
    data[, c(vaVars) := TF$values]
  } else {
    stop("no values given to sum!")
  }
  rm(values)
  
  # additionally, values to compute internal weights by: -----------------------
  iwVar <- NULL
  if (is.character(weights) && 
      pmatch(weights, c("internal", "cohort"), nomatch = 0L)) {
    iw <- substitute(internal.weights.values)
    iw <- evalPopArg(data = origData, iw, DT = TRUE,
                     enclos = PF, recursive = TRUE,
                     types = c("character", "expression", "list", "NULL"))
    
    if (length(iw) > 1L) stop("Argument 'internal.weights.values' ",
                              "must produce only one column.")
    if (length(iw) == 1L && is.character(weights) && 
        pmatch(weights, c("internal", "cohort"), nomatch = 0L)) {
      iwVar <- makeTempVarName(names=c(names(data), names(origData)), pre = "iw_")
      data[, c(iwVar) := TF$iw]
    }
    
    if (length(iwVar) == 0L) {
      stop("Requested computing internal weights, but no values to compute ",
           "internals weights with were supplied (internal error: If you see ",
           "this, complain to the package maintainer).")
    }
    rm(iw)
  }
  
  
  # other category vars to keep ------------------------------------------------
  boSub <- by.other
  by.other <- evalPopArg(data = origData, arg = boSub, DT = TRUE, enclos = enclos)
  boVars <- NULL
  if (length(by.other) > 0) {
    boVars <- names(by.other)
    data[, c(boVars) := TF$by.other]
  } 
  rm(by.other)
  
  # check for aliased columns --------------------------------------------------
  aliased_cols(data, cols = c(prVars, adVars, boVars))
  
  # check for conflicting column names -----------------------------------------
  dupCols <- c(prVars, adVars, boVars, vaVars, iwVar)
  dupCols <- unique(dupCols[duplicated(dupCols)])
  if (length(dupCols) > 0L) {
    dupCols <- paste0("'", dupCols, "'", collapse = ", ")
    stop("Following column names duplicated (columns created by arguments ", 
         "print, adjust, etc.): ", dupCols, ". If you see this, please ensure ",
         "you are not passing e.g. the same column to both for adjusting ",
         "and stratification (printing).")
  }
  
  # check for NA values --------------------------------------------------------
  ## NOTE: NA values of print/by.other are OK. values/adjust are NOT.
  
    NAs <- data[, lapply(.SD, function(x) is.na(x)), .SDcols = c(vaVars, iwVar, adVars)]
    NAs <- rowSums(NAs) > 0L
    if (sum(NAs)) {
      if (!is.null(NA.text)) {
        NA.text <- gsub(x = NA.text, pattern = "%%NA_COUNT%%", 
                        replacement = sum(NAs))
        warning(NA.text)
      }
      data <- data[!NAs]
    }

  # inflate data ---------------------------------------------------------------
  ## on the other hand we aggregate data to levels of print, adjust and 
  ## by.other; on the other hand the data will always have tabulating variables
  ## represented as cross-joined, e.g. merge(1:5, 1:2).
  ## this means some rows might have zeroes as values in the 'values'
  ## columns.
  ## (necessary for correct standardization with weights)
  
  ## NOTE: have to do CJ by hand: some levels of adjust or something may not
  ## have each level of e.g. fot repeated!
  
  sortedLevs <- function(x) {
    if (!is.factor(x)) return(sort(unique(x)))
    
    factor(levels(x), levels(x), levels(x))
  }
  cj <- list()
  cj <- lapply(data[, .SD, .SDcols =  c(prVars, adVars, boVars)], sortedLevs)
  
  ## e.g. if data only has fot = seq(0, 4, 1/12), but want to display table
  ## with fot = seq(0, 5, 1/12). Very important sometimes for handy usage
  ## of weights.
  if (length(custom.levels) > 0) cj[names(custom.levels)] <- custom.levels
  
  
  ## SPECIAL: if e.g. a survival time scale with breaks seq(0, 5, 1/12)
  ## is to be "compressed" to breaks 0:5, and the latter breaks were passed
  ## via custom.levels, the following ensures e.g. intervals between 0 and 1
  ## are aggregated to the same row in what follows after.
  if (!is.null(custom.levels.cut.low)) {
    cl_msg <- paste0("Internal error: tried to cut() variables in  ",
                     "internally used work data that did not exist. ",
                     "If you see this, complain to the ",
                     "package maintainer. Bad variables: %%VARS%%.")
    all_names_present(data, custom.levels.cut.low, msg = cl_msg)
    all_names_present(cj, custom.levels.cut.low, msg = cl_msg)
    
    for (var in custom.levels.cut.low) {
      set(data, j = var, value = cutLow(data[[var]], breaks = cj[[var]]))
    }
    
    ## NOTE: if used cutlow(), then assume passed values via custom.levels
    ## also contained the roof of the values which should not be repeated.
    cj[custom.levels.cut.low] <- lapply(cj[custom.levels.cut.low],
                                        function(elem) elem[-length(elem)])
  }
  
  ## form data.table to merge by - the merge will inflate the data.
  cj <- do.call(function(...) CJ(..., unique = FALSE, sorted = FALSE), cj)
  
  ## inflate & aggregate.
  setkeyv(data, c(prVars, adVars, boVars))
  data <- data[cj, lapply(.SD, sum), .SDcols = c(vaVars, iwVar), by = .EACHI]
  
  for (k in c(vaVars, iwVar)) {
    data[is.na(get(k)), (k) := 0]
  }
  
  setcolsnull(data, tmpDum)
  prVars <- setdiff(prVars, tmpDum); if (length(prVars) == 0) prVars <- NULL
  
  ## merge in weights ----------------------------------------------------------
  
  if (!is.null(weights)) {
    
    if (is.list(weights) && !is.data.frame(weights)) {
      ## in case one of the elements is a standardization scheme string,
      ## such as "world_x_y". 
      whChar <- which(unlist(lapply(weights, is.character)))
      if (sum(whChar)) {
        
        weights[whChar] <- lapply(weights[whChar], function(string) {
          ## expected to return data.frame with 1) age groups 2) weights
          ## as columns.
          stdr.weights(string)[[2]]
        })
        
      }
      
      ## now list only contains numeric weights i hope.
      
    }
    
    ## NOTE: adjust used here to contain levels of adjust arguments only
    adjust <- list()
    if (length(adVars) > 0L) {
      adjust <- lapply(data[, eval(adVars), with = FALSE], sortedLevs)
    }
    
    if (is.character(weights)) {
      
      if (pmatch(weights, c("internal", "cohort"), nomatch = 0L)) {
        
        
        all_names_present(data, iwVar,
                          msg = paste0(
                            "Internal error: expected to have variable ",
                            "%%VARS%% in working data but didn't. Complain ",
                            "to the pkg maintainer if you see this."
                          ))
        
        weights <- mapply(function(levs, colname) {
          setkeyv(data, colname)
          data[.(levs), lapply(.SD, sum), .SDcols = eval(iwVar), by = .EACHI]
        }, 
        colname = names(adjust), levs = adjust,
        SIMPLIFY = FALSE)
        
        weights <- lapply(weights, function(x) {
          x[[iwVar]]
        })
        
        setcolsnull(data, iwVar)
        setkeyv(data, c(prVars, adVars, boVars))
      } else {
        ## expected to return data.frame with 1) age groups 2) weights
        ## as columns.
        weights <- stdr.weights(weights)
        weights <- weights[[2]]
        
      }
    } 
    
    if (!is.data.frame(weights) && is.vector(weights)) {
      ## note: lists are vectors
      if (!is.list(weights)) {
        weights <- list(weights) ## was a vector of values
        setattr(weights, "names", adVars[1])
      }
      
      weVars <- names(weights)
      weights <- weights[names(adjust)]
      
      adjust <- do.call(function(...) CJ(..., unique = FALSE, sorted = FALSE), adjust)
      weights <- do.call(function(...) CJ(..., unique = FALSE, sorted = FALSE), weights)
      
      weVars <- paste0(weVars, ".w")
      setnames(weights, adVars, weVars)
      weights[, (adVars) := adjust]
      
      set(weights, j = "weights", value = 1L)
      for (k in weVars) {
        set(weights, j = "weights", value = weights$weights * weights[[k]])
      }
      setcolsnull(weights, delete = weVars, soft = FALSE)
      
      ## NOTE: weights will be repeated for each level of print,
      ## and for each level of print the weights must sum to one for things
      ## to work.
      weights[, weights := weights/sum(weights)]
      
    }
    
    if (!is.data.frame(weights)) {
      stop("Something went wrong: 'weights' was not collated into a ",
           "data.frame to merge with data. ",
           "Blame the package maintainer please!")
    }
    ## at this points weights is a data.frame.
    weights <- data.table(weights)
    weights[, weights := as.double(weights)]
    
    ## ensure repetition by print levels if some adjust levels
    ## that exist in weights do not exist in data.
    ## NOTE: weights data.frame has at least as many levels as adjust column
    ## in data (or it has more sometimes).
    wm <- lapply(adVars, function(chStr) {
      col <- weights[[chStr]]
      if (is.factor(col)) return(levels(col))
      sort(unique(col))
      })
    names(wm) <- adVars
    if (length(prVars)) {
      wm[prVars] <- lapply(prVars, function(chStr) {
        col <- data[[chStr]]
        if (is.factor(col)) return(levels(col))
        sort(unique(col))
      })
    }
    wm <- do.call(CJ, wm)
    setDT(wm)
    
    weights <- merge(wm, weights, by = adVars, all.x = TRUE, all.y = TRUE)
    
    byCols <- subsetDTorDF(weights, select = prVars)
    if (!length(prVars)) byCols <- NULL
    weights[, weights := weights/sum(weights), by = eval(byCols)]
    rm(byCols)
    
    data <- merge(data, weights, by = c(prVars, adVars), 
                  all.x = TRUE, all.y = FALSE)
    
    if (any(is.na(data$weights))) {
      ## should not be any NAs since we checked for level congruence
      ## in checkWeights
      stop("Internal error: some weights were NA after merging to working ",
           "data. Complain to the package maintainer if you see this.")
    }
    
    
  }
  
  setattr(data, "makeWeightsDT", list(prVars = prVars, adVars = adVars, 
                                      boVars = boVars, vaVars = vaVars, 
                                      NAs = NAs))
  return(data[])
  
}

checkCharWeights <- function(w) {
  if (is.character(w)) {
    if (length(w) != 1L) {
      stop("weights supplied as a character string must be of length one.")
    }
    if (!pmatch(w, c("internal", "cohort"), nomatch = 0L)) {
      stdr.weights(w)
    }
  }
}

checkWeights <- function(weights, adjust) {
  ## INTENTION: given a list/DF/vector/string specifying weights
  ## and a data.frame/list of the adjusting variables,
  ## checks they are congruent and complains if not.
  allowed_classes <- c("list","data.frame","integer","numeric","character",
                       "NULL")
  if (!any(class(weights) %in% allowed_classes)) {
    stop("weights must be either a list, a data.frame, a numeric variable, ",
         "or a character string specifing the weighting scheme to use. ",
         "See ?direct_standardization for more information.") 
  }
  
  if (is.list(weights) && !is.data.frame(weights) && 
      length(adjust) != length(weights)) {
    stop("Mismatch in numbers of variables (NOT necessarily in the numbers of ",
         "levels/values within the variables) in adjust (", length(adjust), 
         " variables) and weights (", length(weights)," variables); ",
         "make sure each given weights vector has a corresponding ",
         "variable in adjust and vice versa. ",
         "See ?direct_standardization for more information.")
  }
  
  if (is.list(weights)) {
    isChar <- unlist(lapply(weights, is.character))
    if (any(isChar)) {
      lapply(weights[isChar], checkCharWeights)
      weights[isChar] <- lapply(weights[isChar], function(string) {
        if (pmatch(string, c("cohort", "internal"), nomatch = 0L)) {
          stop("List of weights had 'cohort' or 'internal' as at least one ",
               "element, which is currently not supported. ",
               "See ?direct_standardization for more information.")
        }
        stdr.weights(string)[[2]]
        })
    }
  }

  if (is.character(weights)) {
    checkCharWeights(weights)
    if (pmatch(weights, c("internal", "cohort"), nomatch = 0L)) {
      ## done checking since internal weights are pretty fool-proof.
      return(invisible())
    }
    ## if not, pass along as vector of weights.
    weights <- stdr.weights(weights)[[2]]
  }
  
  if (is.numeric(weights)) { 
    if (length(adjust) != 1L) {
      stop("Weights is a numeric vector of weights, ",
           "but there are more or less than one adjusting variable. ",
           "See ?direct_standardization for more information.")
    }
    weights <- list(weights)
    names(weights) <- names(adjust)
  }
  
  ## by now either a list or a data.frame of weights...
  adVars <- names(adjust)
  weVars <- names(weights)
  if (is.data.frame(weights)) {
    if (!"weights" %in% weVars) {
      stop("data.frame of weights did not have column named 'weights'. ",
           "see ?direct_standardization for more information.")
    }
    weVars <- setdiff(weVars, "weights")
  }
  
  badAdVars <- setdiff(adVars, weVars)
  badWeVars <- setdiff(weVars, adVars)
  if (length(badAdVars) > 0) {
    stop("Mismatch in names of variables in adjust and weights; ",
         "following adjust variables not mentioned in weights: ", 
         paste0("'", badAdVars, "'", collapse = ", "))
  }
  
  if (length(badWeVars) > 0) {
    stop("Mismatch in names of variables in adjust and weights; ",
         "following weights variables not mentioned in adjust: ",
         paste0("'", badWeVars, "'", collapse = ", "))
  }
  
  if (is.data.frame(weights)) {
    
    levDiff <- lapply(names(adjust), function(var) {
      !all(adjust[[var]] %in% weights[[var]])
    })
    levDiff <- unlist(levDiff)
    if (any(levDiff)) {
      ## take only first conflicting variable for brevity of error message-
      badVar <- names(adjust)[1]
      badLevs <- setdiff(adjust[[badVar]], weights[[badVar]])
      badLevs <- paste0("'", badLevs, "'", collapse = ", ")
      stop("Missing levels in weights data.frame in variable '", badVar, "': ",
           badLevs, ". These levels were found to exist in the corresponding ",
           "adjusting variable. ",
           "Usual suspects: adjusting variable is a factor and you ",
           "only supplied weights for unique values in your data ",
           "as opposed to the levels of the factor, which may contain levels ",
           "that no row has. Try table(yourdata$yourvariable).")
    }
    
  } else {
    weights <- as.list(weights)
    weights <- weights[adVars]
    
    ## check variable levels
    adjust <- lapply(adjust, function(elem) {
      if (is.factor(elem)) {
        levels(elem)
      } else {
        sort(unique(elem))
      }
    })
    
    weLen <- unlist(lapply(weights, length))
    adLen <- unlist(lapply(adjust, length))
    badLen <- names(adjust)[weLen != adLen]
    
    if (length(badLen) > 0) {
      stop("Mismatch in numbers of levels/unique values in adjusting variables ",
           "and lengths of corresponding weights vectors. ",
           "Names of mismatching variables: ", 
           paste0("'", badLen, "'", collapse = ", "), ". There were ",
           weLen[weLen != adLen], " weights and ", adLen[weLen != adLen], 
           " adjusting variable levels.")
    }
  }
  
  
  
  invisible()
}








