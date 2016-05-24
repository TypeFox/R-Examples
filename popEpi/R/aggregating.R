
#' @title Set 'aggre' attributes to an object by modifying in place
#' @author Joonas Miettinen
#' @description Coerces an R object to an \code{aggre} object, identifying
#' the object as one containing aggregated counts, person-years and other
#' information. \code{setaggre} modifies in place without taking any copies.
#' Retains all other attributes.
#' @param x a \code{data.frame} or \code{data.table}
#' @param values a character string vector; the names of value variables
#' @param by a character string vector; the names of variables by which 
#' \code{values} have been tabulated
#' @details 
#' 
#' \code{setaggre} sets \code{x} to the \code{aggre} class in place 
#' without taking a copy as e.g. \code{as.data.frame.XXX} functions do; see e.g. 
#' \code{\link[data.table]{setDT}}.
#' 
#' @seealso 
#' \code{\link{as.aggre}}
#' 
#' @export setaggre
#' @examples 
#' df <- data.frame(sex = rep(c("male", "female"), each = 5), 
#'                  obs = rpois(10, rep(7,5, each=5)), 
#'                  pyrs = rpois(10, lambda = 10000))
#' setaggre(df, values = c("obs", "pyrs"), by = "sex")
setaggre <- function(x, values = NULL, by = setdiff(names(x), values)) {
  ## input: aggregated data in data.frame or data.table format
  ## intention: any user can define their data as an aggregated data set
  ## which will be usable by survtab / sir / other
  ## output: no need to do x <- setaggre(x); instead modifies attributes in place;
  ## sets "aggre.meta" attribute, a list of names of various variables.
  ## survtab for aggregated data will need this attribute to work.
  all_names_present(x, c(values, by))
  
  if (!inherits(x, "aggre")) {
    cl <- class(x)
    wh <- which(cl %in% c("data.table", "data.frame"))
    wh <- min(wh)
    
    ## yes, from zero: in case only one class
    cl <- c(cl[0:(wh-1)], "aggre", cl[wh:length(cl)])
    setattr(x, "class", cl)
  }
  
  
  setattr(x, "aggre.meta", list(values = values, by = by))
  invisible(x)
}

#' @title Coercion to Class \code{aggre}
#' @author Joonas Miettinen
#' @description Coerces an R object to an \code{aggre} object, identifying
#' the object as one containing aggregated counts, person-years and other
#' information. 
#' @param x an R object to coerce to \code{aggre}; must be 
#' a \code{data.frame} or \code{data.table}
#' @param values a character string vector; the names of value variables
#' @param by a character string vector; the names of variables by which 
#' \code{values} have been tabulated
#' @param ... arguments passed to or from methods
#' @seealso 
#' \code{\link{setaggre}} for modifying in place
#' 
#' 
#' @examples 
#' df <- data.frame(sex = rep(c("male", "female"), each = 5), 
#'                  obs = rpois(10, rep(7,5, each=5)), 
#'                  pyrs = rpois(10, lambda = 10000))
#' dt <- as.data.table(df)
#' 
#' df <- as.aggre(df, values = c("pyrs", "obs"), by = "sex")
#' dt <- as.aggre(dt, values = c("pyrs", "obs"), by = "sex")
#' 
#' class(df)
#' class(dt)
#' 
#' @export
as.aggre <- function(x, values = NULL, by = setdiff(names(x), values), ...) {
  UseMethod("as.aggre", x)
}

#' @describeIn as.aggre Coerces a \code{data.frame} to an \code{aggre} object
#' @export
as.aggre.data.frame <- function(x, values = NULL, by = setdiff(names(x), values), ...) {
  x <- copy(x)
  setaggre(x, values = values, by = by, ...)
  setattr(x, "class", c("aggre", "data.frame"))
  x[]
}

#' @describeIn as.aggre Coerces a \code{data.table} to an \code{aggre} object
#' @export
as.aggre.data.table <- function(x, values = NULL, by = setdiff(names(x), values), ...) {
  x <- copy(x)
  setaggre(x, values = values, by = by, ...)
  setattr(x, "class", c("aggre", "data.table", "data.frame"))
  x[]
}

#' @describeIn as.aggre Default method for \code{as.aggre} (stops computations
#' if no class-specific method found)
#' @export
as.aggre.default <- function(x, ...) {
  stop(gettextf("cannot coerce class \"%s\" to 'aggre'", deparse(class(x))), 
       domain = NA)
}



#' @title Aggregation of split \code{Lexis} data
#' @author Joonas Miettinen
#' @description Aggregates a split \code{Lexis} object by given variables 
#' and / or expressions into a long-format table of person-years and 
#' transitions / end-points. Automatic aggregation over time scales
#' by which data has been split if the respective time scales are mentioned
#' in the aggregation argument to e.g. intervals of calendar time, follow-up time
#' and/or age.
#' @param lex a \code{Lexis} object split with e.g. 
#' \code{\link[Epi]{splitLexis}} or \code{\link{splitMulti}}
#' @param by variables to tabulate (aggregate) by.
#' \link[=flexible_argument]{Flexible input}, typically e.g.
#' \code{by = c("V1", "V2")}. See Details and Examples.
#' @param type determines outputted levels to which data is aggregated varying
#' from returning only rows with \code{pyrs > 0} (\code{"unique"}) to
#' returning all possible combinations of variables given in \code{aggre} even
#' if those combinations are not represented in data (\code{"full"}); 
#' see Details
#' @param sum.values optional: additional variables to sum by argument
#'  \code{by}. \link[=flexible_argument]{Flexible input}, typically e.g.
#' \code{sum.values = c("V1", "V2")}
#' @param subset a logical condition to subset by before computations;
#' e.g. \code{subset = area \%in\% c("A", "B")}
#' @param verbose \code{logical}; if \code{TRUE}, the function returns timings
#' and some information useful for debugging along the aggregation process
#' @details 
#' 
#' \strong{Basics}
#' 
#' \code{aggre} is intented for aggregation of split \code{Lexis} data only.
#' See \code{\link[Epi]{Lexis}} for forming \code{Lexis} objects by hand
#' and e.g. \code{\link[Epi]{splitLexis}}, \code{\link{splitLexisDT}}, and
#' \code{\link{splitMulti}} for splitting the data. \code{\link{lexpand}}
#' may be used for simple data sets to do both steps as well as aggregation
#' in the same function call.
#' 
#' Here aggregation refers to computing person-years and the appropriate events
#' (state transitions and end points in status) for the subjects in the data.
#' Hence, it computes e.g. deaths (end-point and state transition) and 
#' censorings (end-point) as well as events in a multi-state setting
#' (state transitions).
#' 
#' The result is a long-format \code{data.frame} or \code{data.table}
#' (depending on \code{options("popEpi.datatable")}; see \code{?popEpi})
#' with the columns \code{pyrs} and the appropriate transitions named as
#' \code{fromXtoY}, e.g. \code{from0to0} and \code{from0to1} depending
#' on the values of \code{lex.Cst} and \code{lex.Xst}.
#' 
#' 
#' \strong{The by argument}
#' 
#' The \code{by} argument determines the length of the table, i.e.
#' the combinations of variables to which data is aggregated.  
#' \code{by} is relatively flexible, as it can be supplied as
#' 
#' \itemize{
#'  \item{a character string vector, e.g. \code{c("sex", "area")}, 
#'  naming variables existing in \code{lex}}
#'  \item{an expression, e.g. \code{factor(sex, 0:1, c("m", "f"))} 
#'  using any variable found in \code{lex}}
#'  \item{a list (fully or partially named) of expressions, e.g. 
#'  \code{list(gender = factor(sex, 0:1, c("m", "f"), area)}}
#' }
#' 
#' Note that expressions effectively allow a variable to be supplied simply as
#' e.g. \code{by = sex} (as a symbol/name in R lingo).
#' 
#' The data is then aggregated to the levels of the given variables 
#' or expression(s). Variables defined to be time scales in the supplied 
#' \code{Lexis} are processed in a special way: If any are mentioned in the
#' \code{by} argument, intervals of them are formed based on the breaks
#' used to split the data: e.g. if \code{age} was split using the breaks 
#' \code{c(0, 50, Inf)}, mentioning \code{age} in \code{by} leads to
#' creating the \code{age} intervals \code{[0, 50)} and \code{[50, Inf)}
#' and aggregating to them. The intervals are identified in the output
#' as the lower bounds of the appropriate intervals.
#' 
#' The order of multiple time scales mentioned in \code{by} matters,
#' as the last mentioned time scale is assumed to be a survival time scale
#' for when computing event counts. E.g. when the data is split by the breaks
#' \code{list(FUT = 0:5, CAL = c(2008,2010))}, time lines cut short at
#' \code{CAL = 2010} are considered to be censored, but time lines cut short at
#' \code{FUT = 5} are not. See Return.
#' 
#' \strong{Aggregation types (styles)}
#' 
#' It is almost always enough to aggregate the data to variable levels
#' that are actually represented in the data 
#' (default \code{aggre = "unique"}; alias \code{"non-empty"}). 
#' For certain uses it may be useful
#' to have also "empty" levels represented (resulting in some rows in output
#' with zero person-years and events); in these cases supplying
#' \code{aggre = "full"} (alias \code{"cartesian"}) causes \code{aggre}
#' to determine the Cartesian product of all the levels of the supplied 
#' \code{by} variables or expressions and aggregate to them. As an example
#' of a Cartesian product, try
#' 
#' \code{merge(1:2, 1:5)}.
#' 
#' @return 
#' A long \code{data.frame} or \code{data.table} of aggregated person-years 
#' (\code{pyrs}), numbers of subjects at risk (\code{at.risk}), and events
#' formatted \code{fromXtoY}, where \code{X} and \code{X} are states 
#' transitioning from and to or states at the end of each \code{lex.id}'s 
#' follow-up (implying \code{X} = \code{Y}). Subjects at risk are computed 
#' in the beginning of an interval defined by any Lexis time scales and 
#' mentioned in \code{by}, but events occur at any point within an interval.
#' 
#' When the data has been split along multiple time scales, the last
#' time scale mentioned in \code{by} is considered to be the survival time 
#' scale with regard to computing events. Time lines cut short by the
#' extrema of non-survival-time-scales are considered to be censored
#' ("transitions" from the current state to the current state).
#' 
#' @seealso \code{\link{aggregate}} for a similar base R solution,
#' and \code{\link{ltable}} for a \code{data.table} based aggregator. Neither
#' are directly applicable to split \code{Lexis} data.
#' 
#' 
#' @examples 
#' 
#' ## form a Lexis object
#' library(Epi)
#' data(sibr)
#' x <- sibr[1:10,]
#' x[1:5,]$sex <- 0 ## pretend some are male
#' x <- Lexis(data = x,
#'            entry = list(AGE = dg_age, CAL = get.yrs(dg_date)),
#'            exit = list(CAL = get.yrs(ex_date)),
#'            entry.status=0, exit.status = status)
#' x <- splitMulti(x, breaks = list(CAL = seq(1993, 2013, 5), 
#'                                  AGE = seq(0, 100, 50)))
#' 
#' ## these produce the same results (with differing ways of determining aggre)
#' a1 <- aggre(x, by = list(gender = factor(sex, 0:1, c("m", "f")), 
#'              agegroup = AGE, period = CAL))
#' 
#' a2 <- aggre(x, by = c("sex", "AGE", "CAL"))
#' 
#' a3 <- aggre(x, by = list(sex, agegroup = AGE, CAL))
#' 
#' ## returning also empty levels
#' a4 <- aggre(x, by = c("sex", "AGE", "CAL"), type = "full")
#' 
#' ## computing also expected numbers of cases
#' x <- lexpand(sibr[1:10,], birth = bi_date, entry = dg_date,
#'              exit = ex_date, status = status %in% 1:2, 
#'              pophaz = popmort, fot = 0:5, age = c(0, 50, 100))
#' x$d.exp <- with(x, lex.dur*pop.haz)
#' ## these produce the same result
#' a5 <- aggre(x, by = c("sex", "age", "fot"), sum.values = list(d.exp))
#' a5 <- aggre(x, by = c("sex", "age", "fot"), sum.values = "d.exp")
#' a5 <- aggre(x, by = c("sex", "age", "fot"), sum.values = d.exp)
#' ## same result here with custom name
#' a5 <- aggre(x, by = c("sex", "age", "fot"), 
#'              sum.values = list(expCases = d.exp))
#'              
#' ## computing pohar-perme weighted figures
#' x$d.exp.pp <- with(x, lex.dur*pop.haz*pp)
#' a6 <- aggre(x, by = c("sex", "age", "fot"), 
#'              sum.values = c("d.exp", "d.exp.pp"))
#' ## or equivalently e.g. sum.values = list(expCases = d.exp, expCases.p = d.exp.pp).
#' @export
aggre <- function(lex, by = NULL, type = c("unique", "full"), sum.values = NULL, subset = NULL, verbose = FALSE) {
  
  allTime <- proc.time()
  
  lex.id <- at.risk <- NULL ## APPEASE R CMD CHECK
  
  PF <- parent.frame(1L)
  TF <- environment()
  
  type <- match.arg(type[1], c("non-empty", "unique", "full", "cartesian"))
  if (type == "cartesian") type <- "full"
  if (type == "non-empty") type <- "unique"
  
  if (verbose) cat("Aggregation type: '", type, "' \n", sep = "")
  
  checkLexisData(lex)
  
  breaks <- copy(attr(lex, "breaks"))
  checkBreaksList(lex, breaks)
  
  allScales <- copy(attr(lex, "time.scales"))
  if (length(allScales) == 0 ) {
    stop("could not determine names of time scales; ",
         "is the data a Lexis object?")
  }
  
  ## subset --------------------------------------------------------------------
  subset <- substitute(subset)
  subset <- evalLogicalSubset(lex, subset)
  
  ## check sum.values ----------------------------------------------------------
  sumSub <- substitute(sum.values)
  sum.values <- evalPopArg(lex[1:min(nrow(lex), 20L), ], arg = sumSub, 
                     enclos = PF, recursive = TRUE, DT = TRUE)
  sumType <- attr(sum.values, "arg.type")
  sumVars <- attr(sum.values, "all.vars")
  sumSub <- attr(sum.values, "quoted.arg")
  if (is.null(sum.values)) {
    sumType <- "NULL"
    sumVars <- NULL
    sumSub <- quote(list())
  }
  badSum <- names(sum.values)[!sapply(sum.values, is.numeric)]
  if (length(badSum) > 0L) {
    badSum <- paste0("'", badSum, "'", collapse = ", ")
    stop("Following variables resulting from evaluating supplied sum.values ",
         "argument are not numeric and cannot be summed: ", badSum, 
         ". Evaluated sum.values: ", deparse(sumSub))
  }
  
  
  ## by argument type -------------------------------------------------------
  ## NOTE: need to eval by AFTER cutting time scales!
  
  ags <- substitute(by)
  if (verbose) cat("Used by argument:", paste0(deparse(ags)),"\n")
  
  ## NOTE: with recursive = TRUE, evalPopArg digs deep enough to find
  ## the actual expression (substituted only once) and returns that and other
  ## things in attributes. Useful if arg substituted multiple times.
  by <- evalPopArg(data = lex[1:min(nrow(lex), 20),], 
                      arg = ags, DT = TRUE, enclos = PF, recursive = TRUE)
  ags <- attr(by, "quoted.arg") 
  av <- attr(by, "all.vars")
  argType <- attr(by, "arg.type")
  
  if (is.null(by)) {
    ags <- substitute(list())
    av <- NULL
    argType <- "NULL"
    type <- "unique"
  }
  if (verbose) cat("Type of by argument:", argType, "\n")
  
  ## take copy of lex ----------------------------------------------------------
  ## if lex is a data.table, this function gets really complicated.
  ## if copy is taken only of necessary vars, it should be fine.
  keepVars <- unique(c("lex.id", allScales, "lex.dur", 
                       "lex.Cst", "lex.Xst", av, sumVars))
  lex.orig <- lex
  lex <- subsetDTorDF(lex, subset = subset, select = keepVars)
  lex <- setDT(copy(lex))
  forceLexisDT(lex, breaks = breaks, allScales = allScales, key = FALSE)
  
  ## ensure no observations outside breaks limits are left in
  lex <- intelliDrop(lex, breaks = breaks)
  
  setkeyv(lex, c("lex.id", allScales[1]))
  setcolsnull(lex, delete = setdiff(allScales, names(breaks)))
  
  ## cut time scales for aggregating if needed ---------------------------------
  aggScales <- intersect(av, allScales)
  if (any(!aggScales %in% names(breaks))) {
    aggScales <- paste0("'", setdiff(aggScales, names(breaks)), "'", collapse = ", ")
    stop("Requested aggregating by time scale(s) by which data ",
         "has not been split: ", aggScales)
  }
  
  ## before cutting, find out which rows count towards "at.risk" figure:
  ## of all scales in aggScales, the last one (or the only one) is assumed
  ## to be the survival time scale.
  tmpAtRisk <- makeTempVarName(lex, pre = "at.risk_")
  set(lex, j = tmpAtRisk, value = TRUE)
  survScale <- NULL
  
  if (length(aggScales) > 0) {
    cutTime <- proc.time()
    ## "at.risk" counts subjects at risk in the beginning of the survival
    ## time scale interval.
    survScale <- aggScales[length(aggScales)]
    lex[, c(tmpAtRisk) := lex[[survScale]] %in% breaks[[survScale]] ]
    catAggScales <- paste0("'", aggScales, "'", collapse = ", ")
    if (verbose) {
      cat("Following time scales mentioned in by argument and will be",
          "categorized into intervals (defined by breaks in object",
          "attributes) for aggregation:", catAggScales, "\n")
    }
    
    ## NEW METHOD: use a copy of lex and just modify in place.
    
    for (sc in aggScales) {
      set(lex, j = sc, value = cutLow(lex[[sc]], breaks = breaks[[sc]]))
    }
    
    if (verbose) cat("Time taken by cut()'ting time scales: ", timetaken(cutTime), "\n")
  }
  
  othVars <- setdiff(av, aggScales)
  if (verbose && length(othVars) > 0) {
    catOthVars <- paste0("'", othVars, "'", collapse = ", ")
    cat("Detected the following non-time-scale variables to be utilized in aggregating:", catOthVars, "\n")
  }
  
  ## eval by -------------------------------------------------------------------
  ## NOTE: needed to eval by AFTER cutting time scales!
  by <- evalPopArg(data = lex, arg = ags, DT = TRUE, enclos = PF, recursive = TRUE)
  byNames <- names(by)
  
  ## computing pyrs ------------------------------------------------------------
  ## final step in determining at.risk:
  ## a lex.id is at.risk only once per by-level
  pyrsTime <- proc.time()
  vdt <- data.table(pyrs = lex$lex.dur, at.risk = lex[[tmpAtRisk]], 
                    lex.id = lex$lex.id)
  pyrs <- vdt[, .(pyrs = sum(pyrs), 
                  at.risk = sum(!duplicated(lex.id) & at.risk)), 
              keyby = by]
  setDT(pyrs)
  
  rm(vdt)
  sumNames <- NULL
  if (sumType != "NULL") {
    if (sumType == "character") {
      sumNames <- evalPopArg(lex, sumSub, n = 1L, DT = FALSE, recursive = TRUE, enclos = PF)
      sum.values <- lex[, lapply(.SD, sum), keyby = by, .SDcols = c(sumNames)]
    } else {
      sum.values <- evalPopArg(lex, sumSub, n = 1L, enclos = PF)
      sumNames <- names(sum.values)
      sumTmpNames <- makeTempVarName(lex, pre = sumNames)
      set(lex, j = sumTmpNames, value = sum.values)
      sum.values <- lex[, lapply(.SD, sum), keyby = by, .SDcols = sumTmpNames]
      setnames(sum.values, sumTmpNames, sumNames)
      setcolsnull(lex, sumTmpNames)
    }
    
    setDT(sum.values)
    pyrs <- merge(pyrs, sum.values, all = TRUE)
    rm(sum.values)
  }
  
  
  if (verbose) cat("Time taken by aggregating pyrs: ", timetaken(pyrsTime), "\n")
  
  valVars <- setdiff(names(pyrs), byNames) ## includes pyrs and anything created by sum
  
  pyrs[is.na(pyrs), pyrs := 0]
  pyrs <- pyrs[pyrs > 0]
  
  aggPyrs <- pyrs[, sum(pyrs)] 
  lexPyrs <- sum(lex.orig$lex.dur[subset])
  pyrsDiff <- aggPyrs - lexPyrs
  if (!isTRUE(all.equal(aggPyrs, lexPyrs, scale = NULL))) {
    warning("Found discrepancy of ", abs(round(pyrsDiff, 4)), " ",
            "in total aggregated pyrs compared to ",
            "sum(lex$lex.dur); compare results by hand and make sure ",
            "settings are right \n")
  }
  rm(subset, aggPyrs, lexPyrs)
  
  ## cartesian output ----------------------------------------------------------
  if (type == "full") {
    carTime <- proc.time()
    
    varsUsingScales <- NULL
    
    ## which variables used one time scale? and which one?
    ## will only be used in cartesian stuff.
    if (argType == "character") {
      varsUsingScales <- intersect(by, aggScales)
      whScaleUsed <- varsUsingScales
    } else if (argType != "NULL") {
      ## note: ags a substitute()'d list at this point always if not char
      whScaleUsed <- lapply(ags[-1], function(x) intersect(all.vars(x), aggScales))
      ## only one time scale should be used in a variable!
      oneScaleTest <- any(sapply(whScaleUsed, function(x) length(x) > 1L))
      if (oneScaleTest) stop("Only one Lexis time scale can be used in any one variable in by argument!")
      varsUsingScales <- byNames[sapply(whScaleUsed, function (x) length(x) == 1L)]
      whScaleUsed <- unlist(whScaleUsed)
    }
    
    ceejay <- lapply(by, function(x) if (is.factor(x)) levels(x) else sort(unique(x)))
    if (length(aggScales) > 0) {
      ## which variables in ceejay used the Lexis time scales from lex?
      
      ceejay[varsUsingScales] <- lapply(breaks[whScaleUsed], function(x) x[-length(x)])
    }
    
    ceejay <- do.call(CJ, ceejay)
    setkeyv(ceejay, byNames)
    setkeyv(pyrs, byNames)
    
    pyrs <- pyrs[ceejay]
    rm(ceejay)
    
    if (verbose) cat("Time taken by making aggregated data large in the cartesian product sense: ", timetaken(carTime), "\n")
  }
  
  
  ## computing events ----------------------------------------------------------
  
  transTime <- proc.time()
  
  if (is.null(by) || (is.data.table(by) && nrow(by) == 0L)) {
    
    by <- quote(list(lex.Cst, lex.Xst))
    
  } else {
    for (var in c("lex.Cst", "lex.Xst")) {
      set(by, j = var, value = lex[[var]])
    }
  }
 
  
  ## NOTE: this will ensure correct detection of censorings:
  ## observations cut short by e.g. period window's edge
  ## will be considered a censoring if the breaks along that time scale
  ## are not passed to detectEvents (assuming the survival time scale is
  ## used in by). If no time scale mentioned in by, then all endings
  ## of observations are either censorings or events.
  detBr <- breaks[survScale]
  if (!length(survScale)) detBr <- NULL
  hasEvent <- detectEvents(lex, breaks = detBr, by = "lex.id") %in% 1:2
  ## is language if user supplied by = NULL 
  if (!is.language(by)) by <- by[hasEvent]
  
  trans <- lex[hasEvent, list(obs = .N), keyby = by]
  
  rm(by, lex)
  
  
  if (verbose) cat("Time taken by aggregating events: ", timetaken(transTime), "\n")
  
  ## casting & merging ---------------------------------------------------------
  
  mergeTime <- proc.time()
  setDT(trans)
  setDT(pyrs)
  
  ## tmpTr to be used in casting
  tmpTr <- makeTempVarName(trans, pre = "trans_")
  trans[, c(tmpTr) := paste0("from", lex.Cst, "to", lex.Xst)]
  transitions <- sort(unique(trans[[tmpTr]]))
  trans[, c("lex.Cst", "lex.Xst") := NULL]
  
  ## note: need tmpDum if by = NULL for correct casting & merging
  tmpDum <- makeTempVarName(trans)
  byNames <- c(byNames, tmpDum)
  byNames <- setdiff(byNames, c("lex.Cst", "lex.Xst"))
  trans[, c(tmpDum) := 1L]
  pyrs[, c(tmpDum) := 1L]
  
  valVars <- unique(c(valVars, transitions))
  
  trans <- cast_simple(trans, rows = byNames, columns = tmpTr, values = "obs")
  
  setkeyv(trans, NULL); setkeyv(pyrs, NULL) ## dcast.data.table seems to keep key but row order may be funky; this avoids a warning
  setkeyv(trans, byNames); setkeyv(pyrs, byNames)
  trans <- trans[pyrs]; rm(pyrs)
  
  trans[, c(tmpDum) := NULL]
  byNames <- setdiff(byNames, tmpDum)
  setcolorder(trans, c(byNames, valVars))
  
  if (verbose) cat("Time taken by merging pyrs & transitions: ", timetaken(mergeTime), "\n")
  
  if (length(valVars) > 0L) {
    trans[, c(valVars) := lapply(.SD, function(x) {
      x[is.na(x)] <- 0
      x
    }), .SDcols = c(valVars)]
  }
  
  
  ## final touch ---------------------------------------------------------------
  setDT(trans)
  alloc.col(trans) ## some problems with internal errors...
  setaggre(trans, values = c("pyrs", "at.risk", transitions, sumNames), by = byNames)
  setattr(trans, "breaks", breaks)
  if (!getOption("popEpi.datatable")) setDFpe(trans)
  if (verbose) cat("Time taken by aggre(): ", timetaken(allTime), "\n")
  
  
  
  trans[]
}












