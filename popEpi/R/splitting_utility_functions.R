all_breaks_in <- function(bl1, bl2, x = NULL) {
  ## INTENTION: return TRUE/FALSE depending on whether bl1 is a subset of bl2;
  ## this means that each element in bl1 exists in bl2, and that those elements
  ## are each subsets of the corresponding elements in bl2.
  ## this is handy to check whether the some Lexis data has already
  ## been split using the breaks in bl1.
  ## NOTE: use checkBreakList() on each list separately before this.
  
  if (!is.list(bl1) || !is.list(bl2)) {
    stop("Arguments bl1 and bl2 must be lists of breaks as supplied to e.g. ",
         "splitMulti.")
  }
  
  if (inherits(x, "Lexis")) {
    checkLexisData(x)
    checkBreaksList(x, bl1)
    checkBreaksList(x, bl2)
  }
  
  ce <- intersect(names(bl1), names(bl2))
  if (length(ce) != length(bl1)) return(FALSE)
  
  test <- mapply(function(l1, l2) {
    all(l1 %in% l2)
  }, l1 = bl1, l2 = bl2[ce], SIMPLIFY = FALSE)
  
  all(unlist(test))
}


checkBreaksList <- function(x, breaks = list(fot = 0:5)) {
  if (is.null(breaks)) stop("breaks is NULL")
  if (!is.list(breaks)) stop("breaks needs to be a list")
  if (!is.data.frame(x)) stop("x needs to be a data.frame")
  timeScales <- names(breaks)
  if (length(breaks) == 0L) stop("length of breaks list is zero")
  if (length(timeScales) != length(breaks)) stop("breaks needs to be a fully named list")
  
  bad_scales <- setdiff(timeScales, names(x))
  if (length(bad_scales) > 0) {
    stop("at least one breaks list name wasn't a variable in data; bad names: ", 
         paste0("'", bad_scales, "'", collapse = ", "))
  }
  lens <- lapply(breaks, function(el) if (is.null(el)) -1 else length(el))
  badLens <- names(lens[unlist(lens) == 0L])
  if (length(badLens)) {
    badLens <- paste0("'", badLens, "'", collapse = ", ")
    stop("Elements in breaks list for the following time scales were of ",
         "length zero but not NULL: ", badLens, ". Breaks list may only ",
         "contain elements of length > 0 or elements that are NULL.")
  }
  invisible(NULL)
}

checkPophaz <- function(lex, ph, haz.name = "haz") {
  ## INTENTION: checks a Lexis data set against the pophaz data set for
  ## consistency (e.g. existing variables to merge by)
  
  if (!is.data.frame(ph)) {
    stop("Data set containing population/expected hazards must be a data.frame",
         " (or a data.table, which is also a data.frame).")
  }
  
  if (!haz.name %in% names(ph)) {
    stop("Data set containing population/expected hazards does not contain a ",
         "column named 'haz'. Make sure the name is exactly that (",
         "case sensitive).")
  }
  
  if (haz.name %in% names(lex)) {
    stop("Lexis data set already contains a column named 'haz', which is a ",
         "reserved name for the population hazard variable to be merged. ",
         "Please rename/delete 'haz' from/in your Lexis data first.")
  }
  
  if (!is.data.frame(ph)) {
    stop("Data set of expected/population hazards must be a data.frame.")
  }
  
  bn <- setdiff(names(ph), haz.name)
  
  if (length(bn) == 0L) {
    stop("No variables in expected/population hazards data set to use in merge ",
         "with Lexis data. Ensure that the pop. haz. data set containts some ",
         "variables to merge by (e.g. sex, calendar year, and age group)")
  }
  if (!all(bn %in% names(lex))) {
    badbn <- paste0("'", setdiff(bn, names(lex)), "'", collapse = ", ")
    stop("Lexis data set did not have following variable(s) that were in ",
         "the expected/population hazards data set: ", badbn,". ",
         "Ensure you have supplied the right data and that the names of the ",
         "intended variables match.")
  }
  
  mergeVars <- setdiff(names(ph), haz.name)
  dup <- any(duplicated(as.data.table(ph), by = mergeVars))
  if (dup) {
    stop("Supplied data set of population/expected hzards has duplicated rows ",
         "by the variables ", paste0("'",mergeVars, "'", collapse = ", "),
         " which prevents correct usage of the data set. Please ensure no rows",
         " area duplicated in the data set before proceeding. Tip: use e.g. ",
         "duplicated(PH, by = c('V1', 'V2')) to check for duplicatedness in ",
         "your data set (here named PH) by the variables V1 and V2."
         )
  }
  
  invisible()
}


globalVariables(c("lex.dur", "lex.Xst", "lex.Cst"))
intelliCrop <- function(x, breaks = list(fot = 0:5), allScales = NULL, cropStatuses = FALSE, tol = .Machine$double.eps^0.5) {
  
  TF <- environment()
  checkBreaksList(x = x, breaks = breaks)
  breaks[unlist(lapply(breaks, length)) == 0L] <- NULL
  if (!is.data.table(x)) stop("x needs to be a data.table")
  
  cropScales <- names(breaks)
  
  all_names_present(x, c("lex.dur", allScales))
  
  if (cropStatuses) {
    ## wrapped in data.table to avoid conflicting variable names with x
    ## (maybe origDur might exist in x)
    origEnd <- x$lex.dur + x[[allScales[1L]]]
  }
  
  
  deltas <- mapply(function(b, y) pmax(min(b), y) - y, SIMPLIFY = FALSE,
                   b = breaks, y = subsetDTorDF(x, select = cropScales))
  ## below: baseline (zero value without assigning zero of bad class)
  deltas <- c(deltas, list(x[[cropScales[1]]][1L] - x[[cropScales[1]]][1L]))
  deltas <- do.call(pmax, deltas)
  
  x[, c(allScales) := .SD + TF$deltas, .SDcols = eval(allScales)]
  x[, lex.dur := lex.dur - TF$deltas]
  
  durs <- mapply(function(b, y) max(b) - y, SIMPLIFY = FALSE,
                 b = breaks, y = subsetDTorDF(x, select = cropScales))
  durs$lex.dur <- x$lex.dur
  durs <- do.call(pmin, durs)
  ## now have max durs by row, i.e. up to roof of breaks at most,
  ## or to ((original lex.dur) - (deltas)) if that is smaller.
  ## (being cropped or exiting before roof of breaks)
  
  x[, lex.dur := TF$durs]
  
  if (cropStatuses) {
    harmonizeStatuses(x, C = "lex.Cst", X = "lex.Xst")
    x[lex.dur + x[[allScales[1L]]] + tol < origEnd, lex.Xst := lex.Cst]
  }
  
  invisible(x)
}

harmonizeStatuses <- function(x, C = "lex.Cst", X = "lex.Xst") {
  
  clC <- class(x[[C]])
  clX <- class(x[[X]])
  tyC <- typeof(x[[C]])
  tyX <- typeof(x[[X]])
  cl <- c(clC, clX)
  
  if (tyC != tyX && clC != clX) {
    if (is.numeric(x[[C]]) && is.numeric(x[[X]])) {
      harmonizeNumeric(x = x, v1="lex.Cst", v2="lex.Xst")
      
    } else if (is.factor(x[[C]]) || is.factor(x[[X]])) {
      if (!is.factor(x[[C]])) set(x, j = C, value = as.factor(x[[C]]))
      if (!is.factor(x[[X]])) set(x, j = X, value = as.factor(x[[X]]))
      
    }
  }
  
  if (any(cl == "factor")) {
    harmonizeFactors(x = x,  v1="lex.Cst", v2="lex.Xst")
  }
  
}

harmonizeNumericTimeScales <- function(x, times = NULL) {
  ## INTENTION: given a Lexis data set with some time scales, ensure
  ## that the classes of the time scales comply to the lowest denominator,
  ## e.g. "double" and "integer" -> both "double"
  
  if (is.null(times)) {
    times <- c(attr(x, "time.scales"), "lex.dur")
  }
  
  msg <- paste0("Expected working data to have time scales %%VARS%%, but it ",
                "didn't. This is an internal error: If you see this, complain ",
                "to the package maintainer.")
  all_names_present(x, times, msg = msg)
  xt <- lapply(times, function(ch) x[[ch]])
  names(xt) <- times
  
  harmoClasses <- c("numeric", "integer", "difftime")
  cl <- lapply(xt, class)
  wh <- unlist(lapply(cl, function(ch) {
    any(ch %in% harmoClasses)
  }))
  ha <- times[wh]
  hacl <- unique(unlist(cl[wh]))
  
  if (length(ha) > 1L) {
    ## more than one class present and need to use common lowest denom
    newMode <- as.double
    
    if (all(ha %in% c("integer", "difftime"))) {
      ## all numeric times are integers or difftimes
      newMode <- as.integer
    }
    for (var in ha) {
      ## modify in place
      set(x, j = var, value = newMode(x[[var]]))
    }
    
    
  }
  invisible(NULL)
}

harmonizeNumeric <- function(x, v1="lex.Cst", v2="lex.Xst") {
  ## assumes v1, v2 are numeric variable names in x  
  
  if (!is.numeric(x[[v1]]) || !is.numeric(x[[v2]])) {
    print(class(x[[v1]]))
    print(class(x[[v2]]))
    stop("v1 and/or v2 is/are not of class numeric")
  }
  
  if (!is.integer(x[[v1]])) set(x, j = v1, value = try2int(x[[v1]]))
  if (!is.integer(x[[v2]])) set(x, j = v2, value = try2int(x[[v2]]))
  
  if (typeof(x[[v1]]) != typeof(x[[v2]])) {
    
    if (is.double(x[[v1]]))  set(x, j = v1, value = as.double(x[[v1]]))
    if (is.double(x[[v2]]))  set(x, j = v2, value = as.double(x[[v2]]))
    
  }  
  
}

harmonizeFactors <- function(x, v1="lex.Cst", v2="lex.Xst") {
  ## assumes v1, v2 are factor names in x
  
  if (!is.factor(x[[v1]]) || !is.factor(x[[v2]])) {
    stop("v1 and/or v2 is/are not of class factor")
  }
  
  glab1 <- union(levels(x[[v1]]), levels(x[[v2]]))
  glab2 <- union(levels(x[[v2]]), levels(x[[v1]]))
  
  
  
  setattr(x[[v1]], "levels", glab1)
  setattr(x[[v2]], "levels", glab2)
  
}

intelliDrop <- function(x, breaks = list(fot = 0:5), dropNegDur = TRUE, check = FALSE, tol = .Machine$double.eps^0.5, subset = NULL)  {
  
  if (!is.data.table(x)) {
    stop("x needs to be a data.table; if you see this message, complain ",
         "to the package maintainer")
  }
  checkBreaksList(x = x, breaks = breaks)
  breaks[unlist(lapply(breaks, length)) == 0L] <- NULL
  timeScales <- names(breaks)
  
  if (check) {
    checkLexisData(x)
  }
  
  ra <- lapply(breaks, range)
  ra <- lapply(ra, diff)
  ts <- names(sort(unlist(ra))) ## shortest first
  mi <- lapply(breaks, min)
  ma <- lapply(breaks, max)
  
  substi <- substitute(subset)
  subset <- evalLogicalSubset(x, substiset = substi)
  
  if (dropNegDur) subset[subset] <- subset[subset] & x$lex.dur[subset] > 0L
  
  e <- environment()
  
  ## figure out latest exit and first entry; don't need to test for dropping
  ## if e.g. all left follow-up before the max in breaks
  max_end <- lapply(ts, function(ch) min(x[[ch]] + x$lex.dur))
  min_start <- lapply(ts, function(ch) max(x[[ch]]))
  names(max_end) <- names(min_start) <- ts
  
  for (k in ts) {
    mik <- mi[[k]]
    mak <- ma[[k]]
    
    if (max_end[[k]] < mak + tol) {
      tmpSD <- x[e$subset, .SD, .SDcols = c(e$k, "lex.dur")]
      tmpSD <- setDT(lapply(tmpSD, as.numeric))
      subset[subset] <- rowSums(tmpSD)  <=  e$mak + e$tol
    }
    if (min_start[[k]] + tol > mik) {
      subset[subset] <- x[e$subset, .SD[[e$k]]  > e$mik - e$tol, 
                          .SDcols = c(e$k)]
    }
    
    if (all(!subset)) {
      stop("Dropped all remaining rows from data when subsetting by the  ",
           "Lexis time scale '", k, "'. Range of values in data: ", 
           paste0(round(range(x[[k]]),4), collapse = "-"), ". Min/Max breaks ",
           "(used to subset data): ", mik, "/", mak, ".")
    }
    
  }
  
  
  x[e$subset, ]
}


matchBreakTypes <- function(lex, breaks, timeScale, modify.lex = FALSE) {
  if (is.character(breaks)) {
    breaks <- as.IDate(breaks)
  }
  clb <- class(breaks)
  clb <- clb[length(clb)]
  cts <- class(lex[[timeScale]])
  cts <- cts[length(cts)]
  
  if (clb != cts) {
    if (is.Date(breaks) && !is.Date(lex[[timeScale]])) {
      breaks <- try2int(as.double(breaks))
    } else if (is.integer(breaks) && is.double(lex[[timeScale]])) {
      breaks <- as.double(breaks)
    } else if (is.double(breaks) && is.integer(lex[[timeScale]])) {
      breaks <- try2int(breaks)
    }
    
  }
  
  if (modify.lex && clb != cts) {
    if (!is.Date(breaks) && is.Date(lex[[timeScale]])) {
      
      if (clb == "double") {
        set(lex, j = timeScale, value = as.double(lex[[timeScale]]))
      } else {
        set(lex, j = timeScale, value = as.integer(lex[[timeScale]]))
      }
      
    } else if (is.double(breaks) && is.integer(lex[[timeScale]])) {
      set(lex, j = timeScale, value = as.double(lex[[timeScale]]))
    }
    
  }
  breaks
}

protectFromDrop <- function(breaks, lower = FALSE) {
  old_breaks <- copy(breaks)
  if (length(breaks) == 0L) {
    stop("Length of breaks to 'protect' from dropping is zero.")
  }
  if (is.Date(breaks)) {
    breaks <- c(breaks, max(breaks) + 1e4L)
    if (lower) breaks <- c(min(breaks) - 1e4L, breaks)
    
  } else if (is.integer(breaks))  {
    breaks <- c(breaks, 1e6L)
    if (lower) breaks <- c(-1e6L, breaks)
    
  } else if (is.double(breaks)) {
    breaks <- c(breaks, Inf)
    if (lower) breaks <- c(-Inf, breaks)
    
  } else {
    stop("breaks were not Date, integer or double")
  }
  setattr(breaks, "unprotected", old_breaks)
  breaks
}

unprotectFromDrop <- function(breaks) {
  up <- attr(breaks, "unprotected")
  if (is.null(up) || length(up) == 0L) {
    stop("Could not 'unprotect' breaks from dropping as the required ",
         "attribute was not found. If you see this it is most likely ",
         "an internal error and you should complain to the pkg maintainer.")
  }
  up
}




setLexisDT <- function(data, entry, exit, entry.status, exit.status, id = NULL, select = NULL) {
  
  if (!is.data.table(data)) stop("not a data.table")
  if (inherits(data, "Lexis")) stop("already a Lexis object")
  
  if (!is.null(select) && !is.character(select)) stop("select was not a character vector of names")
  
  entry <- substitute(entry)
  exit <- substitute(exit)
  entry <- eval(entry, envir = data, enclos = parent.frame())
  exit <- eval(exit, envir = data, enclos = parent.frame())
  enNames <- names(entry)
  exNames <- names(exit)
  
  timeScales <- union(enNames, exNames)
  if (any(timeScales %in% names(data))) stop("at least one named time scales already present in data; original names mandatory")
  enNeeded <- setdiff(timeScales, enNames)
  enPresent <- setdiff(timeScales, enNeeded)
  durVar <- intersect(enNames, exNames)
  if (length(durVar) > 1) stop("you have more than 1 time scales in both entry and exit; only one mandatory")
  
  enVars <- paste0(enNames, "_en")
  exVars <- paste0(exNames, "_ex")
  setattr(entry, "names", enVars)
  setattr(exit, "names", exVars)
  
  l <- as.data.table(c(entry, exit))
  rm(entry, exit)
  
  ## duration
  exV <- paste0(durVar, "_ex")
  enV <- paste0(durVar, "_en")
  set(l, j = "lex.dur", value = l[[exV]] - l[[enV]])
  rm(exV, enV)
  
  ## time scale starting points
  if (length(enNeeded) > 0) {
    for (ts in enNeeded) {
      exV <- paste0(ts, "_ex")
      set(l, j = ts, value = l[[exV]] - l$lex.dur)
    }
  }
  setnames(l, paste0(enPresent, "_en"), enPresent)
  
  # no longer need time scale end points
  for (k in exVars) {
    set(l, j = k, value = NULL)
  }
  
  ## status definition
  data[, lex.Cst := entry.status]
  data[, lex.Xst := exit.status]
  
  harmonizeStatuses(data, C = "lex.Cst", X = "lex.Xst")
  
  ## all time scales etc. into data
  data[, names(l) := l]
  
  
  id <- substitute(id)
  id <- eval(id, envir = data, enclos = parent.frame())
  if (!is.null(id)) set(data, j = "lex.id", value = id)
  rm(id)
  
  if (!is.null(select)) {
    
    delVars <- setdiff(names(data), c(names(l), select))
    if (length(delVars) > 0) {
      l[, (delVars) := NULL]
    }
  }
  
  rm(l)
  lexVars <- c("lex.id", timeScales, "lex.dur", "lex.Cst", "lex.Xst")
  setcolorder(data, c(lexVars, setdiff(names(data), lexVars)))
  
  setattr(data, "time.scales", timeScales)
  setattr(data, "time.since", rep("", times = length(timeScales)))
  setattr(data, "class", c("Lexis", "data.table", "data.frame"))
  
  
}

checkLexisData <- function(lex, check.breaks = FALSE) {
  ## INTENTION: checks Lexis attributes
  ## OUTPUT: nothing
  
  if (is.null(lex) || nrow(lex) == 0) stop("Data is NULL or has zero rows")
  if (!inherits(lex, "Lexis")) stop("Data not a Lexis object")
  allScales <- attr(lex, "time.scales")
  if (length(allScales) == 0) stop("no time scales appear to be defined; is data a Lexis object?")
  
  badScales <- setdiff(allScales, names(lex))
  if (length(badScales) > 0) {
    badScales <- paste0("'", badScales, "'", collapse = ", ")
    stop("Following time scales found in data's attributes but not present in data: ", badScales)
  }
  
  lexVars <- c("lex.dur", "lex.id", "lex.Cst", "lex.Xst")
  blv <- setdiff(lexVars, names(lex))
  if (length(blv) > 0) {
    blv <- paste0("'", blv, "'", collapse = ", ")
    stop("Following Lexis variables not found in data: ", blv)
  }
  
  if (check.breaks) {
    BL <- attr(lex, "breaks")
    if (is.null(BL)) stop("No breaks list in data attributes")
    checkBreaksList(lex, breaks = BL)
  }
  
  invisible()
}


splitMultiPreCheck <- function(data = NULL, breaks = NULL, ...) {
  
  ## INTENTION: checks for discrepancies between data and breaks, etc.
  ## OUTPUT: cleaned-up list of breaks
  checkLexisData(data)
  allScales <- attr(data, "time.scales")
  
  if (!is.null(breaks) && !is.list(breaks)) stop("breaks must be a list; see examples in ?splitMulti")
  if (is.null(breaks)) {
    breaks <- list(...)
    breaks <- breaks[intersect(names(breaks), allScales)]
  }
  
  if (length(breaks) == 0) stop("no breaks defined!")
  
  splitScales <- names(breaks)
  ## NULL breaks imply not used
  for (k in splitScales) {
    if (length(breaks[[k]]) == 0) {
      breaks[k] <- NULL
    }
  } 
  
  checkBreaksList(x = data, breaks = breaks)
  
  splitScales <- names(breaks)
  
  if (!all(splitScales %in% allScales)) {
    stop("breaks must be a list with at least one named vector corresponding to used time scales \n
         e.g. breaks = list(fot = 0:5)")
  }
  
  if (!all(splitScales %in% names(data))) {
    stop("At least one vector name in breaks list is not a variable name in the data")
  }
  breaks
}

forceLexisDT <- function(x, breaks = NULL, allScales = NULL, key = TRUE) {
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  setattr(x, "breaks", breaks)
  setattr(x, "time.scales", allScales)
  # alloc.col(x)
  if (key) setkeyv(x, c("lex.id", names(breaks)[1L]))
  invisible(x)
}


doCutLexisDT <- function(lex, cut = dg_date, timeScale = "per", by = "lex.id", n = 1L) {
  
  checkLexisData(lex, check.breaks = FALSE)
  
  x <- unique(lex, by = by)
  cut <- evalq(cut, envir = lex, enclos = parent.frame(n = n + 1L))
  delta <- cut - x[[timeScale]]
  
  allScales <- attr(lex, "time.scales")
  
  setDT(x)
  for (v in allScales) {
    set(x, j = v, value = x[[v]] + delta)
  }
  
  set(x, j = "lex.dur", value = 0)
  
  tmp <- list()
  tmp$isCut <- makeTempVarName(lex, pre = "isCut_")
  
  set(x, j = tmp$isCut, value = 1L)
  on.exit(setcolsnull(x, unlist(tmp)))
  
  x <- rbindlist(list(lex, x), use.names = TRUE, fill = TRUE)
  x[1:nrow(lex), (tmp$isCut) := 0L]
  
  ## NOTE: new cut row being the first or last row
  ## implies it resides outside old observations
  ## OR it is equal to lowest/highest value
  setkeyv(x, c(by, allScales))
  setkeyv(x, by)
  x <- x[!((duplicated(x) | duplicated(x, fromLast = TRUE)) & get(tmp$isCut) == 0L)]
  stop("not ready")
}

# data <- data.table(birth = 2000:2000, entry=2002:2003, 
#                    exit=2011:2012, event=c(2010,2011), 
#                    status=1:0)
# 
# lex <- lexpand(data = data, birth = birth, entry = entry, 
#                exit = exit, event = event, 
#                id = 1L, entry.status = 99L,
#                status = status, overlapping = T)

lexpile <- function(lex, by = "lex.id", subset = NULL) {
  ## PURPOSE: given several rows per id in a Lexis object,
  ## collate data into form where
  ## - no subject has any overlapping time lines
  ## - lex.Cst and lex.Xst are logical, i.e. 0 -> 1, 1 -> 1, 1 -> 2
  ## this should be made to work with both split and unsplit Lexis data.
  
  data <- NULL # R CMD CHECK appeasement
  
  checkLexisData(lex, check.breaks = FALSE)
  
  allScales <- attr(lex, "time.scales")
  sc <- allScales[1L]

  all_names_present(lex, by)
  
  if (is.character(lex$lex.Cst) || is.character(lex$lex.Xst)) {
    stop("This function requires lex.Cst and lex.Xst to be integer, double (i.e. numeric) or factor variables to determine the order of possible statuses!")
  }
  
  ## need to take copy eventually ----------------------------------------------
  attrs <- attributes(lex)
  subset <- evalLogicalSubset(data, substitute(subset))
  x <- lex[subset,]
  forceLexisDT(x, breaks = attrs$breaks, allScales = attrs$time.scales)
  alloc.col(x)
  
  
  ## ensure status harmony -----------------------------------------------------
  harmonizeStatuses(x = x, X = "lex.Xst", C = "lex.Cst")
  exStat <- if (is.factor(x$lex.Xst)) levels(x$lex.Xst) else sort(unique(x$lex.Xst))
  enStat <- if (is.factor(x$lex.Cst)) levels(x$lex.Cst) else sort(unique(x$lex.Cst))
  allStat <- c(setdiff(enStat, exStat), exStat) ## enStat & exStat equal if factors used
  
  ## avoiding side effects -----------------------------------------------------
  oldKey <- key(lex)
  tmp <- list()
  tmp$order<- makeTempVarName(x, pre = "order_")
  
  on.exit({
    if (length(oldKey) > 0) setkeyv(x, oldKey) else 
      setorderv(x, tmp$order)
  }, add = TRUE)
  
  on.exit({
    setcolsnull(x, unlist(tmp$order), soft = TRUE)
  }, add = TRUE)
  
  x[, c(tmp$order) := 1:.N]
  
  ## check for need for lexpiling ----------------------------------------------
  setkeyv(x, by)
  if (sum(duplicated(x)) == 0L) return(lex)
  
  
  ## figure out what statuses are used -----------------------------------------
  
  tmp$ev <- makeTempVarName(x, pre = "event_")
  x[, c(tmp$ev) := detectEvents(x, breaks = attrs$breaks, by = by)]
  
  tmp$scEnds <- paste0(allScales, "_end")
  tmp$scEnds <- makeTempVarName(lex, pre = tmp$scEnds)
  x[, c(tmp$scEnds) := lapply(.SD, function(x) x + lex$lex.dur), .SDcols = allScales]
  
  ## NOTE: rows for a given subject ending in simultaneously with at least
  ## one being a transition will not be allowed.
  setkeyv(x, c(by, tmp$scEnds[1L]))
  
  whDup <- duplicated(x, fromLast = FALSE) | duplicated(x, fromLast = TRUE)
  dupTest <- x[whDup, 1L %in% unique(.SD), .SDcols = tmp$ev]
  rm(whDup)
  
  if (dupTest) stop("At least one subject had at least two simultaneous events.")
  ## NOTE: if interval ends AND status are the very same for M rows,
  ## then the M rows are necessarily nested with one or more covering
  ## the whole time line. Only need to keep the one.
  setkeyv(x, c(tmp$scEnds, "lex.Cst", "lex.Xst"))
  setorderv(x, c(allScales,tmp$scEnds, "lex.Cst", "lex.Xst"))
  
  x <- unique(x)
  stop("unfinished")
  
}

contractLexis <- function(x, breaks, drop = TRUE) {
  stop("This doesnt do anything yet")
  ## INTENTION: given a Lexis object and breaks,
  ## ensures data is split by the breaks and contracts the split rows
  ## so that the data is split at the level of the supplied breaks.
  ## e.g. with x split by fot = seq(0, 5, 1/12) and with supplying
  ## breaks = list(fot = 0:5), rows within 0-1 are collated into one row etc.
  
  ## PROBLEM: a subject may have e.g. rows spaning 0.5 - 1 which are requested
  ## to be contracted to one row spanning 0-1.
  
  
}




#' @title Prepare Exposure Data for Aggregation
#' @description \code{prepExpo} uses a \code{Lexis} object of periods of exposure
#' to fill gaps between the periods and overall entry and exit times without
#' accumulating exposure time in periods of no exposure, and splits the
#' result if requested.
#' @param lex a \code{\link[Epi]{Lexis}} object with ONLY periods of exposure
#' as rows; one or multiple rows per subject allowed
#' @param freezeScales a character vector naming \code{Lexis} time scales of exposure
#' which should be frozen in periods where no exposure occurs (in the gap
#' time periods) 
#' @param cutScale the \code{Lexis} time scale along which the subject-specific
#' ultimate entry and exit times are specified
#' @param entry an expression; the time of entry to follow-up which may be earlier, at, or after
#' the first time of exposure in \code{freezeScales}; evaluated separately
#' for each unique combination of \code{by}, so e.g. with 
#' \code{entry = min(Var1)} and \code{by = "lex.id"} it 
#' sets the \code{lex.id}-specific minima of \code{Var1} to be the original times
#' of entry for each \code{lex.id}
#' @param exit the same as \code{entry} but for the ultimate exit time per unique
#' combination of \code{by}
#' @param by a character vector indicating variable names in \code{lex},
#' the unique combinations of which identify separate subjects for which
#' to fill gaps in the records from \code{entry} to \code{exit};
#' for novices of \code{{\link{data.table}}}, this is passed to a 
#' \code{data.table}'s \code{by} argument.
#' @param breaks a named list of breaks; 
#' e.g. \code{list(work = 0:20,per = 1995:2015)}; passed on to 
#' \code{\link{splitMulti}} so see that function's help for more details
#' @param freezeDummy a character string; specifies the name for a dummy variable
#' that this function will create and add to output which 
#' identifies rows where the \code{freezeScales} are frozen and where not
#' (\code{0} implies not frozen, \code{1} implies frozen);
#' if \code{NULL}, no dummy is created
#' @param subset a logical condition to subset data by before computations;
#' e.g. \code{subset = sex == "male"}
#' @param verbose logical; if \code{TRUE}, the function is chatty and returns
#' some messages and timings during its run.
#' @param ... additional arguments passed on to \code{\link{splitMulti}}
#' @details 
#' 
#' \code{prepExpo} is a convenience function for the purpose of eventually aggregating 
#' person-time and events in categories of not only normally progressing 
#' \code{Lexis} time scales but also some time scales which should not
#' progress sometimes. For example a person may work at a production facility
#' only intermittently, meaning exposure time (to work-related substances 
#' for example) should not progress outside of periods of work. This allows for
#' e.g. a correct aggregation of person-time and events by categories of cumulative
#' time of exposure.
#' 
#' Given a \code{Lexis} object containing rows (time lines)
#' where a subject is exposed to something (and NO periods without exposure),
#' fills any gaps between exposure periods for each unique combination of \code{by}
#' and the subject-specific "ultimate" \code{entry} and \code{exit} times,
#' "freezes" the cumulative exposure times in periods of no exposure,
#' and splits data using \code{breaks} passed to \code{\link{splitMulti}}
#' if requested. Results in a (split) \code{Lexis} object where \code{freezeScales}
#' do not progress in time periods where no exposure was recorded in \code{lex}.
#' 
#' This function assumes that \code{entry} and \code{exit} arguments are the
#' same for each row within a unique combination of variables named in \code{by}.
#' E.g. with \code{by = "lex.id"} only each \code{lex.id} has a unique value
#' for \code{entry} and \code{exit} at most.
#' 
#' The supplied \code{breaks} split the data using \code{splitMulti}, with
#' the exception that breaks supplied concerning any frozen time scales
#' ONLY split the rows where the time scales are not frozen. E.g.
#' with \code{freezeScales = "work"}, 
#' \code{breaks = list(work = 0:10, cal = 1995:2010)} splits all rows over
#' \code{"cal"} but only non-frozen rows over \code{"work"}.
#' 
#' Only supports frozen time scales that advance and freeze contemporaneously:
#' e.g. it would not currently be possible to take into account the cumulative
#' time working at a facility and the cumulative time doing a single task
#' at the facility, if the two are not exactly the same. On the other hand
#' one might use the same time scale for different exposure types, supply them
#' as separate rows, and identify the different exposures using a dummy variable.
#' @return 
#' 
#' Returns a \code{Lexis} object that has been split if \code{breaks} is specified.
#' The resulting time is also a \code{data.table} if 
#' \code{options("popEpi.datatable") == TRUE} (see: \code{?popEpi})
#' 
#' @import data.table
#' @export
prepExpo <- function(lex, freezeScales = "work", cutScale = "per", entry = min(get(cutScale)),
                     exit = max(get(cutScale)), by = "lex.id", breaks = NULL, freezeDummy = NULL, subset = NULL,
                     verbose = FALSE, ...) {
  
  if (verbose) allTime <- proc.time()
  
  ## check breaks & data -------------------------------------------------------
  breaks <- evalq(breaks)
  dumBreaks <- structure(list(c(-Inf, Inf)), names = cutScale, internal_prepExpo_dummy = TRUE)
  if (is.null(breaks)) breaks <- dumBreaks
  breaks <- splitMultiPreCheck(data = lex, breaks = breaks)
  if (!is.null(attr(breaks, "internal_prepExpo_dummy"))) breaks <- NULL
  checkLexisData(lex)
  oldBreaks <- attr(lex, "breaks")
  if (!is.null(breaks)) checkBreaksList(lex, breaks)
  checkBreaksList(lex, oldBreaks)
  
  
  ## data ----------------------------------------------------------------------
  
  subset <- evalLogicalSubset(data = lex, substitute(subset))
  x <- if (!all(subset)) evalq(lex)[subset, ] else copy(evalq(lex))
  
  setDT(x)
  
  allScales <- attr(lex, "time.scales")
  linkScales <- setdiff(allScales, freezeScales)
  othScales <- setdiff(linkScales, cutScale)
  
  setkeyv(x, c(by, cutScale))
  
  l <- list() ## will hold temp var names; this avoids collisions with names of vars in x
  l$cutScale <- cutScale
  l$freezeScales <- freezeScales
  l$liquidScales <- setdiff(allScales, freezeScales)
  l$by <- by
  rm(cutScale, freezeScales, by)
  
  ## args ----------------------------------------------------------------------
  if (verbose) argTime <- proc.time()
  tol <- .Machine$double.eps^0.75
  
  if (is.character(freezeDummy) && freezeDummy %in% names(lex)) stop("Variable named in freezeDummy already exists in data; freezeDummy is inteded for creating a new dummy for identifying the rows where freezeScales are frozen. Please supply an original variable name to freezeDummy")
  
  if (!is.character(l$by)) stop("by must be given as a vector of character strings naming columns in lex")
  all_names_present(lex, l$by)
  
  enSub <- substitute(entry)
  exSub <- substitute(exit)
  
  PF <- parent.frame(1L)
  l$en <- makeTempVarName(x, pre = "entry_")
  l$ex <- makeTempVarName(x, pre = "exit_")
  x[, c(l$ex, l$en) := list(eval(exSub, envir = .SD, enclos = PF), 
                                eval(enSub, envir = .SD, enclos = PF)), by = c(l$by)]
  
  ## tests disabled for now...
#   testTime <- proc.time()
#   
#   test <- x[, .N, by = list(r = get(l$cutScale) + lex.dur > get(l$ex) - tol)]
#   if(test[r == TRUE, .N] > 0) stop("exit must currently be higher than or equal to the maximum of cutScale (on subject basis defined using by); you may use breaks instead to limit the data")
#   
#   test <- x[, .N, by = list(r = get(l$cutScale) + tol < get(l$en))]
#   if(test[r == TRUE, .N] > 0) stop("entry must currently be lower than or equal to the minimum of cutScale (on subject basis defined using by); you may use breaks instead to limit the data")
#   if (verbose) cat("Finished checking entry and exit. Time taken: ", timetaken(argTime), "\n")
  if (verbose) cat("Finished evaluating entry and exit and checking args. Time taken: ", timetaken(argTime), "\n")
  
  ## create rows to fill gaps --------------------------------------------------
  if (verbose) fillTime <- proc.time()
  x2 <- copy(x)
  x2[, (l$freezeScales) := NA]
  x2 <- rbind(x2, unique(x2, by = c(l$by), fromLast = TRUE))
  
  l$delta <- makeTempVarName(x2, pre = "delta_")
  x2[, (l$delta) := c(get(l$en)[1], get(l$cutScale)[-c(1,.N)], max(get(l$cutScale)+lex.dur)) - get(l$cutScale), by = c(l$by)]
  x2[, c(linkScales) := lapply(mget(linkScales), function(x) x + get(l$delta)), by = c(l$by)]
  
  setcolsnull(x2, l$delta)
  
  l$order <- makeTempVarName(x, pre = "order_")
  x[, (l$order) := (1:.N)*2, by = c(l$by)]
  x2[, (l$order) := (1:.N)*2-1, by = c(l$by)]
  
  x <- rbindlist(list(x, x2))
  rm(x2)
  setkeyv(x, c(l$by, l$cutScale, l$order))
  setkeyv(x, c(l$by))
  set(x, j = l$order, value = as.integer(x[[l$order]]))
  x[, (l$order) := 1:.N, by = c(l$by)]
  
  if (verbose) cat("Finished expanding data to accommodate filling gaps. Time taken: ", timetaken(fillTime), "\n")
  ## handle time scale values --------------------------------------------------
  if (verbose) valueTime <- proc.time()
  
  l$CSE <- makeTempVarName(x, pre = paste0(l$cutScale, "_end_"))
  l$LCS <- makeTempVarName(x, pre = paste0("lead1_",l$cutScale, "_"))
  x[, (l$CSE) := lex.dur + get(l$cutScale)]
  x[, (l$LCS)  := shift(get(l$cutScale), n = 1L, type = c("lead"), fill = NA), by = c(l$by)]
  
  
  x[!duplicated(x, fromLast = TRUE), c(l$LCS, l$CSE) := get(l$ex)]
  x[, (l$CSE) := pmin(get(l$LCS), get(l$CSE))]
  x[, (l$cutScale) := sort(c(get(l$en)[1L],shift(get(l$CSE), n = 1L, type = "lag", fill = NA)[-1])), by = c(l$by)]
  x[, lex.dur := get(l$CSE) - get(l$cutScale)]
  
  ## bring up other than frozen and cut scales to bear -------------------------
  x[, (othScales) := lapply(mget(othScales), function(x) {min(x) + c(0, cumsum(lex.dur)[-.N])}), by = c(l$by)]
  
  
  ## frozen scales should make sense cumulatively ------------------------------
  ## indicates frozenness: 0 = not frozen, 1 = frozen
  l$frz <- makeTempVarName(x, pre = "frozen_")
  x[, (l$frz) := 0L]
  frozens <- x[,is.na(get(l$freezeScales[1]))]
  x[frozens, (l$frz) := 1L]
  
  
  ## alternate method: just use lex.durs and only cumulate in non-frozen rows
  x[, (l$freezeScales) := lapply(mget(l$freezeScales), function(x) {
    x <- max(0, min(x-lex.dur, na.rm=TRUE))
    x <- x + c(0, as.double(cumsum(as.integer(!get(l$frz))*lex.dur))[-.N])
  }), by = c(l$by)]
  
  x <- x[lex.dur > .Machine$double.eps^0.5, ]
  
  if (verbose) cat("Finished computing correct values for time scales. Time taken: ", timetaken(valueTime), "\n")
  
  ## splitting separately ------------------------------------------------------
  if (!is.null(breaks)) {
    if (verbose) splitTime <- proc.time()
    x_frozen <- x[get(l$frz) == 1L,]
    x <- x[get(l$frz) == 0L]
    forceLexisDT(x, allScales = allScales, breaks = oldBreaks)
    
    ## NOTE: since we only split by the frozen time scales by pretending
    ## they are NOT Lexis time scales (temporarily), and since one should 
    ## pass the appropriate breaks info of pre-existing breaks to splitMulti,
    ## choose only breaks for non-frozen time scales to include in x_frozen's
    ## attributes here.
    ## (e.g. when work history is no longer accumulating)
    frzBreaks <- breaks[l$liquidScales]
    oldFrzBreaks <- oldBreaks[l$liquidScales]
    emptyFrzBreaks <- vector("list", length = length(l$liquidScales))
    names(emptyFrzBreaks) <- l$liquidScales
    if (length(oldFrzBreaks)) {
      emptyFrzBreaks[names(oldFrzBreaks)] <- oldFrzBreaks
    }
    forceLexisDT(x_frozen, allScales = l$liquidScales, breaks = emptyFrzBreaks)
    
    if (length(frzBreaks) > 0) {
      ## do (also) split for all time scales where also the frozen
      ## time scales are split. This is allowed for times where the
      ## frozen time scales have not been frozen
      ## (e.g. work history is accumulating)
      x_frozen <- splitMulti(x_frozen, breaks = frzBreaks, ...)
    }
    
    ## do (also) split where also split
    x <- splitMulti(x, breaks = breaks, ...)
    breaks <- attr(x, "breaks") ## new breaks appended by splitMulti
    
    setDT(x)
    setDT(x_frozen)
    x <- rbindlist(list(x, x_frozen), use.names = TRUE); rm(x_frozen)
    forceLexisDT(x, breaks = breaks, allScales = allScales)
    if (verbose) cat("Finished splitting data. Time taken: ", timetaken(splitTime), "\n")
  }
  
  ## final touch ---------------------------------------------------------------
  
  setDT(x)
  if (is.character(freezeDummy)) setnames(x, l$frz, freezeDummy)
  setkeyv(x, c(l$by, l$order))
  delCols <- setdiff(names(l), c("by", "cutScale", "freezeScales", 
                                 "liquidScales",
                                 "linkScales", "allScales", "othScales"))
  delCols <- unlist(l[delCols])
  setcolsnull(x, keep = names(lex), colorder = TRUE)
  
  setattr(x, "time.scales", allScales)
  setattr(x, "breaks", breaks)
  setattr(x, "time.since", rep("", length(allScales)))
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  if (getOption("popEpi.datatable") == FALSE) setDFpe(x)
  
  if (verbose) cat("Finished prepExpo run. Time taken: ", timetaken(allTime), "\n")
  
  x[]
}



doComparisonWithEpi <- function(lexDT, lexDTdrop, lexDF, breaks) {
  BL <- NULL
  if (!is.list(breaks)) stop("breaks needs to be a list")
  requireNamespace("Epi")
  requireNamespace("testthat")
  
  allScales <- attr(lexDF, "time.scales")
  sc1 <- allScales[1]
  setDT(lexDT)
  setDT(lexDTdrop)
  setDT(lexDF)
  setkeyv(lexDT, c("lex.id", sc1))
  setkeyv(lexDTdrop, c("lex.id", sc1))
  setkeyv(lexDF, c("lex.id", sc1))
  
  testthat::expect_equal(attr(lexDT, "time.scales"), attr(lexDF, "time.scales"))
  testthat::expect_equal(attr(lexDT, "time.since"), attr(lexDF, "time.since"))
  
  testthat::expect_equal(attr(lexDTdrop, "time.scales"), attr(lexDF, "time.scales"))
  testthat::expect_equal(attr(lexDTdrop, "time.since"), attr(lexDF, "time.since"))
  
  doTestBarrage(dt1 = lexDT, dt2 = lexDF, allScales = allScales)
  rm(lexDT)
  
  lexDF <- intelliDrop(x = lexDF, breaks = breaks)
  
  doTestBarrage(dt1 = lexDTdrop, dt2 = lexDF, allScales = allScales)
  
}

doTestBarrage <- function(dt1, dt2, allScales, testTimes = TRUE, testStatuses = TRUE) {
  requireNamespace("Epi")
  requireNamespace("testthat")
  
  lex.id <- NULL ## APPEASE R CMD CHECK
  
  testthat::expect_equal(sum(dt1$lex.dur), 
                         sum(dt2$lex.dur), 
                         check.attributes = FALSE)
  testthat::expect_equal(dt1[, sum(lex.dur), keyby = lex.id]$V1, 
                         dt2[, sum(lex.dur), keyby = lex.id]$V1, 
                         check.attributes = FALSE)
  
  all_names_present(dt1, allScales)
  all_names_present(dt2, allScales)
  
  if (testTimes) {
    for (k in allScales) {
      testthat::expect_equal(dt1[[k]], dt2[[k]], 
                             check.attributes = TRUE)
    }
  }
  
  if (testStatuses) {
    testthat::expect_equal(dt1$lex.Cst, dt2$lex.Cst, check.attributes = FALSE)
    testthat::expect_equal(dt1$lex.Xst, dt2$lex.Xst, check.attributes = FALSE)
    
    testthat::expect_equal(levels(dt1$lex.Cst), levels(dt2$lex.Cst), check.attributes = FALSE)
    testthat::expect_equal(levels(dt1$lex.Xst), levels(dt2$lex.Xst), check.attributes = FALSE)
    
    testthat::expect_true(all(class(dt2$lex.Cst) %in% class(dt1$lex.Cst)))
    testthat::expect_true(all(class(dt2$lex.Xst) %in% class(dt1$lex.Xst)))
  }
  
  invisible(NULL)
}

compareSLDTWithEpi <- function(data, breaks, timeScale) {
  requireNamespace("Epi")
  requireNamespace("testthat")
  
  if (!inherits(data, "Lexis")) stop("data gotta be a Lexis object broseph")
  
  lexDT <- splitLexisDT(data, breaks = breaks, timeScale = timeScale, merge = TRUE, drop = FALSE)
  lexDTdrop <- splitLexisDT(data, breaks = breaks, timeScale = timeScale, merge = TRUE, drop = TRUE)
  lexDF <- splitLexis(data, breaks = breaks, time.scale = timeScale) ## without dropping
  ## this treatment done in splitLexisDT (difftime -> integer -> double)
  harmonizeNumericTimeScales(lexDF, times = c(Epi::timeScales(lexDF), "lex.dur"))
  
  BL <- list(breaks)
  setattr(BL, "names", timeScale)
  
  doComparisonWithEpi(lexDT = lexDT, lexDTdrop = lexDTdrop, lexDF = lexDF, breaks = BL)
  
  invisible(NULL)
}

splitMultiEpi <- function(data, breaks = list(fot = 0:5), drop) {
  
  for (k in names(breaks)) {
    data <- splitLexis(data, breaks = breaks[[k]], time.scale = k)
  }
  
  if (drop) data <- intelliDrop(data, breaks = breaks)
  setDT(data)
  
  data
}

compareSMWithEpi <- function(data, breaks = list(fot=0:5)) {
  requireNamespace("Epi")
  requireNamespace("testthat")
  
  lexDT <- splitMulti(data, breaks = breaks, merge = TRUE, drop = FALSE)
  lexDTdrop <- splitMulti(data, breaks = breaks, merge = TRUE, drop = TRUE)
  lexDF <- splitMultiEpi(data, breaks = breaks, drop = FALSE)
  
  doComparisonWithEpi(lexDT=lexDT, lexDTdrop = lexDTdrop, lexDF=lexDF, breaks = breaks)
  
  invisible(NULL)
}
