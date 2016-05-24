#' @title Split case-level observations
#' @author Joonas Miettinen
#' @description Split a \code{Lexis} object along multiple time scales
#' with speed and ease
#' @param data a Lexis object with event cases as rows
#' @param breaks a list of named numeric vectors of breaks; see Details and Examples
#' @param ... alternate way of supplying breaks as named vectors;
#' e.g. \code{fot = 0:5} instead of \code{breaks = list(fot = 0:5)};
#' if \code{breaks} is not \code{NULL}, \code{breaks} is used and any breaks
#' passed through \code{...} are NOT used
#' @param drop logical; if \code{TRUE}, drops all resulting rows 
#' after expansion that reside outside the time window
#' defined by the given breaks
#' @param merge logical; if \code{TRUE}, retains all variables 
#' from the original data - i.e. original variables are
#' repeated for all the rows by original subject
#' @param verbose logical; if \code{TRUE}, the function is chatty 
#' and returns some messages along the way
#' 
#' 
#' @details 
#' 
#' \code{splitMulti} is in essence a \pkg{data.table} version of
#'  \code{splitLexis} or \code{survSplit} for splitting along multiple
#'  time scales.
#' It requires a Lexis object as input.
#' 
#' The \code{breaks} must be a list of named vectors of the appropriate type. 
#' The breaks are fully explicit and
#' left-inclusive and right exclusive, e.g. \code{fot=c(0,5)} 
#' forces the data to only include time between
#' \code{[0,5)} for each original row (unless \code{drop = FALSE}). 
#' Use \code{Inf} or \code{-Inf} for open-ended intervals,
#'  e.g. \code{per=c(1990,1995,Inf)} creates the intervals 
#'  \code{[1990,1995), [1995, Inf)}.
#'  
#' Instead of specifying \code{breaks}, one may make use of the \code{...}
#' argument to pass breaks: e.g. 
#' 
#' \code{splitMulti(x, breaks = list(fot = 0:5))} 
#' 
#' is equivalent to
#' 
#' \code{splitMulti(x, fot = 0:5)}.
#' 
#' Multiple breaks can be supplied in the same manner. However, if both
#' \code{breaks} and \code{...} are used, only the breaks in \code{breaks}
#' are utilized within the function. 
#' 
#' The \code{Lexis} time scale variables can be of any arbitrary 
#' format, e.g. \code{Date},
#' fractional years (see \code{\link[Epi]{cal.yr}}) and \code{\link{get.yrs}},
#' or other. However, using \code{date} variables (from package \pkg{date})
#' are not recommended, as \code{date} variables are always stored as integers,
#' whereas \code{Date} variables (see \code{?as.Date}) are typically stored
#' in double ("numeric") format. This allows for breaking days into fractions
#' as well, when using e.g. hypothetical years of 365.25 days.
#'  
#' @return
#' A \code{data.table} or \code{data.frame} 
#' (depending on \code{options("popEpi.datatable")}; see \code{?popEpi}) 
#' object expanded to accommodate split observations.
#' 
#' @examples
#' #### let's prepare data for computing period method survivals
#' #### in case there are problems with dates, we first 
#' #### convert to fractional years.
#' \dontrun{
#' library(Epi)
#' x <- Lexis(data=sire, entry = list(fot=0, per=get.yrs(dg_date), age=dg_age), 
#'            exit=list(per=get.yrs(ex_date)), exit.status=status)
#' x2 <- splitMulti(x, breaks = list(fot=seq(0, 5, by = 3/12), per=c(2008, 2013)))
#' # equivalently:
#' x2 <- splitMulti(x, fot=seq(0, 5, by = 3/12), per=c(2008, 2013))
#' 
#' ## using dates; note: breaks must be expressed as dates or days!
#' x <- Lexis(data=sire, entry = list(fot=0, per=dg_date, age=dg_date-bi_date), 
#'            exit=list(per=ex_date), exit.status=status)
#' BL <- list(fot = seq(0, 5, by = 3/12)*365.242199,
#'            per = as.Date(paste0(c(1980:2014),"-01-01")),
#'            age = c(0,45,85,Inf)*365.242199)
#' x2 <- splitMulti(x, breaks = BL, verbose=TRUE)
#' 
#' ## multistate (healty - sick - dead)
#' ## pretend some observation never got cancer
#' set.seed(1L)
#' 
#' sire2 <- copy(sire)
#' sire2$status <- factor(sire2$status, levels = 0:2)
#' levels(sire2$status) <- c("healthy", "dead", "dead")
#' 
#' not_sick <- sample.int(nrow(sire2), 6000L, replace = FALSE)
#' sire2[not_sick, ]$dg_date <- NA
#' sire2[!is.na(dg_date) & status == "healthy", ]$status <- "sick"
#' 
#' xm <- Lexis(data=sire2, entry = list(fot=0, per=get.yrs(bi_date), age=0), 
#'             exit=list(per=get.yrs(ex_date)), exit.status=status)
#' xm2 <- cutLexis(xm, cut = get.yrs(xm$dg_date), timescale = "per", new.state = "sick")
#' xm2[xm2$lex.id == 6L, ]
#' 
#' xm2 <- splitMulti(xm2, breaks = list(fot = seq(0,150,25)))
#' xm2[xm2$lex.id == 6L, ]
#' }
#' 
#' @import data.table 
#' @import Epi
#' 
#' @export splitMulti
#' 
#' @seealso
#' \code{\link[Epi]{splitLexis}}, \code{\link[Epi]{Lexis}},  
#' \code{\link[survival]{survSplit}}, \code{\link{splitLexisDT}}
#' 
splitMulti <- function(data,
                       breaks = NULL,
                       ...,
                       drop=TRUE,
                       merge=TRUE,
                       verbose=FALSE) {
  
  lex.id <- NULL ## APPEASE R CMD CHECK
  
  ## basic checks --------------------------------------------------------------
  if (verbose) {stime <- proc.time()}
  
  breaks <- splitMultiPreCheck(data = data, breaks = breaks, ...)
  
  
  ## check if even need to do splitting ----------------------------------------
  allScales <- attr(data, "time.scales")
  splitScales <- names(breaks)
  
  oldBreaks <- attr(data, "breaks")
  tryCatch(checkBreaksList(data, oldBreaks), error = function(e) {
    stop("Error in splitMulti: \n",
         "Old breaks existing in Lexis data did not pass testing. Error ",
         "message from test: \n", e, call. = FALSE)
  })
  
  ## only do split if all breaks are NOT in the breaks that the data
  ## has already been split by.
  do_split <- TRUE
  do_split <- !all_breaks_in(breaks, oldBreaks, x = data)
  
  if (!do_split) {
    l <- setDT(copy(data))
    setkeyv(l, c("lex.id", allScales[1]))
  } else {
    
    ## prepare data  -------------------------------------------------------------
    
    tmpID <- makeTempVarName(data = data, pre = "TEMP_SPLITMULTI_ID_")
    IDT <- data.table(lex.id = data$lex.id, temp.id = 1:nrow(data), key = "temp.id")
    
    set(data, j = "lex.id", value = 1:nrow(data))
    on.exit(set(data, j = "lex.id", value = IDT[, ]$lex.id))
    
    l <- vector(mode = "list", length = length(splitScales))
    setattr(l, "names", splitScales)
    for (v in splitScales) {
      l[[v]] <- splitLexisDT(data, breaks = breaks[[v]], 
                             merge = TRUE, drop = FALSE, timeScale = v)
      breaks[[v]] <- attr(l[[v]], "breaks")[[v]]
    }
    l <- rbindlist(l)
    on.exit()
    set(data, j = "lex.id", value = IDT$lex.id)
    
    ## lex.id is here 1:nrow(data)
    setkey(l, lex.id)
    setnames(l, "lex.id",  tmpID)
    set(l, j = "lex.id", value = IDT[.(l[[tmpID]]), ]$lex.id)
    rm(data, IDT)
    
    
    if (!merge) setcolsnull(l, keep = c("lex.id", "lex.dur", allScales, "lex.Cst", "lex.Xst", tmpID), soft = FALSE)
    
    v1 <- splitScales[1]
    
    if (length(splitScales) > 1L) {
      
      tmpIE <- makeTempVarName(data = l, pre = "TEMP_SPLITMULTI_intEnd")
      tmpLD <- makeTempVarName(data = l, pre = "TEMP_SPLITMULTI_lagDur")
      l[, (tmpIE) := get(v1) + lex.dur]
      setkeyv(l, c(tmpID, tmpIE))
      
      l <- unique(l)
      ## time scale minima and lex.dur as cumulative duration --------------------
      l[, (allScales) := lapply(.SD, min), .SDcols = allScales, by = c(tmpID)]
      l[, lex.dur := get(tmpIE) - get(v1)]
      l[, (tmpLD) := c(0, lex.dur[-.N]), by = c(tmpID)]
      
      for (k in allScales) {
        set(l, j = k, value = l[[k]] + l[[tmpLD]])
      }
      # non-cumulative duration
      set(l, j = "lex.dur", value = l$lex.dur - l[[tmpLD]])
      set(l, j = tmpLD, value = NULL)
      set(l, j = tmpIE, value = NULL)
    }
    
    setkeyv(l, c(tmpID, v1))
    set(l, j = tmpID, value = NULL)
    
  }
  
  l <- l[lex.dur > 0]
  if (drop) l <- intelliDrop(l, breaks = breaks, dropNegDur = FALSE)
  
  if (nrow(l) == 0) {
    warning("no data left after dropping; check breaks?")
  }
  
  order <- c("lex.id", "lex.multi", allScales, "lex.dur", "lex.Cst", "lex.Xst")
  order <- c(order, setdiff(names(l), order))
  order <- intersect(order, names(l))
  setcolorder(l, order)
  
  if (verbose) cat("time taken by splitting process: ", timetaken(stime), "\n")
  
  setattr(l, "time.scales", allScales)
  setattr(l, "time.since", rep("", times=length(allScales)))
  setattr(l, "breaks", breaks)
  setattr(l, "class", c("Lexis","data.table","data.frame"))
  if (!getOption("popEpi.datatable")) setDFpe(l)
  
  l[]
  
}

globalVariables(".")
