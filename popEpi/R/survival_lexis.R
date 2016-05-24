# library(survival)
# dt[, fot := as.integer(ex_date-dg_date)/365.24]
# s <- dt[, Surv(time = rep(0, nrow(dt)), time2 = fot, event = status %in% 1:2)]

#' @describeIn survtab_ag survtab
#' @export
#' @importFrom survival Surv
survtab <- function(formula, data, adjust = NULL, breaks = NULL, 
                    pophaz = NULL, weights = NULL, surv.type = "surv.rel", 
                    surv.method = "hazard", relsurv.method = "e2", 
                    subset = NULL, 
                    conf.level = 0.95, 
                    conf.type = "log-log",
                    verbose = FALSE) {
  
  TF <- environment()
  PF <- parent.frame()
  this_call <- match.call()
  startTime <- proc.time()
  
  if (!requireNamespace("survival")) {
    stop("Need package 'survival' to proceed")
  }
  
  
  ## checks --------------------------------------------------------------------
  
  if (missing(formula)) stop("Formula not defined!")
  
  checkLexisData(data)
  
  allScales <- attr(data, "time.scales")
  splitScales <- names(breaks)
  
  ## ensure breaks make sense --------------------------------------------------
  oldBreaks <- attr(data, "breaks")
  checkBreaksList(data, breaks = oldBreaks)
  testOldBreaks <- setdiff(oldBreaks, list(NULL))
  if (is.null(breaks) && !length(testOldBreaks)) {
    stop("No breaks supplied via argument 'breaks', and data has not been ",
         "split in advance. Please supply a list of breaks ",
         "to argument 'breaks'")
  }
  if (is.null(breaks)) breaks <- oldBreaks
  checkBreaksList(data, breaks = breaks)
  ## match break types to time scale types
  ## (don't try to match time scales to breaks)
  splitScales <- names(breaks)
  for (k in splitScales) {
    breaks[[k]] <- matchBreakTypes(data, breaks = breaks[[k]], timeScale = k)
  }
  
  comp_pp <- FALSE
  drop <- TRUE
  if (surv.type == "surv.rel" && relsurv.method == "pp") comp_pp <- TRUE
  if (comp_pp) drop <- FALSE
  
  
  ## data & subset -------------------------------------------------------------
  subset <- evalLogicalSubset(data, substitute(subset))
  x <- data[subset, ]; rm(subset)
  setDT(x)
  forceLexisDT(x, breaks = NULL, allScales = allScales, key = TRUE)
  
  ## pre-eval of print & adjust ------------------------------------------------
  
  adSub <- substitute(adjust)
  adTest <- evalRecursive(adSub, env = x, enc = PF)
  if (!is.null(adTest)) {
    adSub <- adTest$argSub
    adVars <- all.vars(adSub)
  } else {
    adSub <- substitute(NULL)
    adVars <- NULL
  } 
  
  formula <- evalRecursive(formula, env = TF, enc = PF)$arg
  foVars <- all.vars(formula)
  
  if (!inherits(formula,"formula")) {
    stop("Argument 'formula' is not a formula object. Usage: e.g. ",
         "Surv(fot, lex.Xst %in% 1:2) ~ sex")
  }
  if (length(formula) != 3L) {
    stop("Argument 'formula'must be two-sided. Usage: e.g. ",
         "Surv(fot, lex.Xst %in% 1:2) ~ sex")
  }
  ## eval print & adjust -------------------------------------------------------
  ## this adjust passed to resulting data's attributes at the end
  adSub <- substitute(adjust)
  adjust <- evalPopArg(data = x, arg = adSub, 
                       enclos = PF, DT = TRUE, 
                       recursive = TRUE)
  
  l <- usePopFormula(form = formula, adjust = adjust, data = x, enclos = PF, 
                     Surv.response = "either")
  prVars <- names(l$print)
  adVars <- names(l$adjust)
  
  
  ## check weights makes sense with respect to adjust --------------------------
  if (length(adVars) > 0L && !is.null(weights)) {
    checkWeights(weights, adjust = l$adjust)
    
  }
  
  ## check pophaz --------------------------------------------------------------
  
  if (surv.type %in% c("surv.rel", "cif.rel")) {
    checkPophaz(x, pophaz, haz.name = "haz")
  }
  pophazVars <- setdiff(names(pophaz), "haz")
  
  ## only keep necessary variables ---------------------------------------------
  
  setcolsnull(x, keep = c("lex.id", "lex.dur", allScales, 
                          "lex.Cst", "lex.Xst", pophazVars))
  if (length(prVars)) x[, c(prVars)] <- l$print
  if (length(adVars)) x[, c(adVars)] <- l$adjust
  
  
  ## simplify event and censoring indicators -----------------------------------
  cens.values <- event.values <- NULL
  all.values <- if (is.factor(l$y$status)) levels(l$y$status) else 
    sort(unique(l$y$status))
  cens.values <- all.values[1L]
  event.values <- setdiff(all.values, cens.values)
  
  if (is.numeric(l$y$status) && all(unique(l$y$status) %in% 0:1)) {
    ## this should apply to situations where status coded 0/1
    ## and both 0/1 present or only 1 present
    if (all(unique(l$y$status) %in% 0L)) {
      stop("All status values were zero, i.e. all obs were censored. ",
           "Check that you passed the correct status variable or --- if this ",
           "was intended --- code the status variable to 0/1 so that 1 ",
           "corresponds to the event taking place and 0 not.")
    }
    cens.values <- 0L
    event.values <- 1L
  }
  
  x[, lex.Cst := NULL]
  x[, lex.Cst := TF$cens.values]
  x[, lex.Xst := NULL]
  x[, lex.Xst := l$y$status]
  harmonizeStatuses(x, C = "lex.Cst", X = "lex.Xst")
  
  if (!surv.type %in% c("cif.obs", "surv.cause")) {
    ## this simplifies computations
    
    x[, lex.Cst := NULL]
    x[, lex.Cst := 0L]
    setcolorder(x, c(intersect(names(data), names(x)), 
                     setdiff(names(x), names(data))))
    
    x[, lex.Xst := as.integer(lex.Xst %in% TF$event.values)]
    cens.values <- 0L
    event.values <- 1L
    if (x[, sum(lex.Xst)] == 0L) {
      stop("There are no events in the data. Ensure that the event argument ",
           "used in Surv() makes sense.")
    }
  }
  
  ## detect which time scale used ----------------------------------------------
  
  survScale <- detectSurvivalTimeScale(lex = x, values = l$y$time)
  
  
  ## crop data to speed up computations ----------------------------------------
  cropBreaks <- breaks
  if (surv.type == "surv.rel" && relsurv.method == "pp")  {
    ## pp-weights have to be computed from entry to follow-up till roof of breaks;
    ## can only crop along the survival time scale
    cropBreaks <- breaks[1L]
    cb <- protectFromDrop(cropBreaks[[1L]], lower = TRUE)
    cb <- c(min(cb), max(cropBreaks[[1L]]))
    cropBreaks[[1L]] <- cb
  }
  
  
  intelliCrop(x = x, breaks = cropBreaks, allScales = allScales, cropStatuses = TRUE)
  x <- intelliDrop(x, breaks = cropBreaks, dropNegDur = TRUE, check = TRUE)
  setDT(x)
  forceLexisDT(x, breaks = oldBreaks, allScales = allScales, key = TRUE)
  
  ## splitting -----------------------------------------------------------------
  
  splitTime <- proc.time()
  setDT(x)
  forceLexisDT(x, breaks = oldBreaks, allScales = allScales, key = TRUE)
  x <- splitMulti(x, breaks = breaks, drop = FALSE, merge = TRUE)
  setDT(x)
  forceLexisDT(x, breaks = breaks, allScales = allScales, key = TRUE)
  if (verbose) cat("Time taken by splitting Lexis data: ", timetaken(splitTime), "\n")
  
  ## pophaz merge --------------------------------------------------------------
  if (!is.null(pophaz)) {
    hazTime <- proc.time()
    haz <- NULL ## appease R CMD CHECK
    x <- cutLowMerge(x, pophaz, by = pophazVars, 
                     mid.scales = intersect(pophazVars, allScales))
    setDT(x)
    forceLexisDT(x, breaks = breaks, allScales =allScales, key = TRUE)
    if (verbose) cat("Time taken by merging population hazards with split Lexis data: ", timetaken(hazTime), "\n")
  }
  
  ## pp computation ------------------------------------------------------------
  ppNames <- d.pp <- d.pp.2 <- d.exp.pp <- ptime.pp <- 
    at.risk.pp <- n.cens.pp <- NULL
  
  if (comp_pp) {
    ppTime <- proc.time()
    setkeyv(x, c("lex.id", survScale))
    comp_pp_weights(x, surv.scale = survScale, 
                    breaks = breaks[[survScale]], haz = "haz", 
                    style = "delta", verbose = FALSE)
    setDT(x)
    forceLexisDT(x, breaks = breaks, allScales = allScales, key = TRUE)
    if (verbose) cat("Time taken by computing Pohar-Perme weights: ", timetaken(ppTime), "\n")
    
    intelliCrop(x = x, breaks = breaks, allScales = allScales, cropStatuses = TRUE)
    x <- intelliDrop(x, breaks = breaks, dropNegDur = TRUE, check = TRUE)
    forceLexisDT(x, breaks = breaks, allScales = allScales, key = TRUE)
    
    ppTime <- proc.time()
    pp <- comp_pp_weighted_figures(x, haz = "haz", pp = "pp", by = "lex.id")
    ppNames <- makeTempVarName(x, pre = names(pp))
    x[, c(TF$ppNames) := TF$pp] ## note: TF$pp avoids conflicts
    rm(pp)
    
    d.pp.2 <- ppNames[substr(ppNames, 1, 13) == "from0to1.pp.2"]
    d.pp <- ppNames[substr(ppNames, 1, 11) == "from0to1.pp"]
    d.pp <- setdiff(d.pp, d.pp.2)
    d.exp.pp <- ppNames[substr(ppNames, 1, 8) == "d.exp.pp"]
    ptime.pp <- ppNames[substr(ppNames, 1, 8) == "ptime.pp"]
    n.cens.pp <- ppNames[substr(ppNames, 1, 11) == "from0to0.pp"]
    n.cens.pp <- n.cens.pp[substr(n.cens.pp, 1,13) != "from0to0.pp.2"]
    at.risk.pp <-  ppNames[substr(ppNames, 1, 10) == "at.risk.pp"]
    d.exp.pp <-  ppNames[substr(ppNames, 1, 8) == "d.exp.pp"]
    
    if (verbose) cat("Time taken by computing Pohar-Perme weighted ",
                     "counts and person-times: ", timetaken(ppTime), "\n")
  }
  
  d.exp <- NULL
  if (surv.type %in% c("surv.rel", "cif.rel") && "haz" %in% names(x)) {
    d.exp <- makeTempVarName(x, pre = "d.exp_")
    x[, c(TF$d.exp) := lex.dur * haz]
  }
  
  ## aggregation ---------------------------------------------------------------
  aggreTime <- proc.time()
  
  ## this includes time scale to compute survivals over
  aggreVars <- c(prVars, adVars, survScale) 
  
  setDT(x)
  forceLexisDT(x, breaks = breaks, allScales = allScales, key = TRUE)
  if (verbose) cat("** verbose messages from aggre(): \n")
  x <- aggre(x, by = aggreVars, verbose = verbose,
             sum.values = c(d.exp, ppNames))
  if (verbose) cat("** end of  verbose messages from aggre() \n")
  setDT(x)
  setattr(x, "class", c("aggre", "data.table", "data.frame"))
  if (verbose) cat("Time taken by aggregating split Lexis data: ", 
                   timetaken(aggreTime), "\n")
  
  ## neater column names -------------------------------------------------------
  ## in case there are zero obs that are censored
  censCols <- paste0("from", cens.values, "to", cens.values)
  if (all(!censCols %in% names(x))) {
    x[, c(censCols) := 0L]
  }
  
  ## e.g. fromAlivetoDead -> Dead; looks better in survtab_ag output
  evCols <- paste0("from", cens.values, "to", c(cens.values, event.values))
  whEC <- which(evCols %in% names(x))
  
  if (sum(whEC)) {
    setnames(x, evCols[whEC], 
             as.character(c(cens.values, event.values)[whEC]))
  }
  
  ## survtab_ag ----------------------------------------------------------------
  dn <- intersect(event.values, names(x))
  if (length(dn) == 0L) {
    stop("Internal error: no event variables in work data. Complain to the ",
         "package maintainer if you see this - unless there are no events ",
         "in the data?")
  }
  n.cens <- intersect(cens.values, names(x))
  
  if (length(prVars) == 0L) {
    prVars <- "1"
  } 
  
  form <- as.formula(paste0(survScale, " ~ ", paste0(prVars, collapse = " + ")))
  
  if (verbose) cat("** verbose messages from survtab_ag(): \n")
  
  st <- survtab_ag(data = x, 
                   formula = TF$form,
                   adjust = TF$adVars,
                   
                   weights = TF$weights, 
                   
                   d = TF$dn, pyrs = "pyrs", n = "at.risk",
                   d.exp = TF$d.exp, n.cens = TF$n.cens,
                   
                   n.pp = TF$at.risk.pp,
                   d.pp = TF$d.pp, d.exp.pp = TF$d.exp.pp, d.pp.2 = TF$d.pp.2, 
                   n.cens.pp = TF$n.cens.pp, pyrs.pp = TF$ptime.pp,
                   
                   surv.type = surv.type,
                   surv.method = surv.method,
                   relsurv.method = relsurv.method,
                   
                   conf.type = conf.type,
                   conf.level = conf.level,
                   
                   verbose = verbose)
  if (verbose) cat("** end of verbose messages from survtab_ag() \n")
  ## attributes ----------------------------------------------------------------
  
  attributes(st)$survtab.meta$call <- this_call
  attributes(st)$survtab.meta$arguments$adjust <- adjust
  attributes(st)$survtab.meta$arguments$conf.type <- conf.type
  attributes(st)$survtab.meta$arguments$conf.level <- conf.level
  
  attributes(st)$survtab.meta$arguments$surv.type <- surv.type
  attributes(st)$survtab.meta$arguments$surv.method <- surv.method
  attributes(st)$survtab.meta$arguments$relsurv.method <- relsurv.method
  
  if (verbose) cat("Total time taken by survtab: ", timetaken(startTime), "\n")
  st
}
# library(Epi)
# library(popEpi)
# dt <- copy(sire)[dg_date < ex_date,]
# dt[, agegr := cut(dg_age, c(0,50,75,Inf))]
# dt[, sex := rbinom(n = .N, size = 1, prob = 0.5)]
# dt <- Lexis(data = dt, entry = list(FUT = 0, AGE = dg_age, CAL = get.yrs(dg_date)),
#             exit = list(CAL = get.yrs(ex_date)), entry.status = 0L, exit.status = status, merge = TRUE)
# pm <- copy(popEpi::popmort)
# setnames(pm, c("agegroup", "year"), c("AGE", "CAL"))
# st <- survtab(data = dt, formula = Surv(FUT, lex.Xst) ~ 1, #adjust = "agegr", 
#                   # pophaz = pm,
#                   surv.type = "surv.obs",
#                   # weights = list(agegr = c(0.2,0.4,0.4)),
#                   breaks = list(FUT = seq(0,5,1/12)))
# st <- survtab(dt, print = NULL, #adjust = "agegr", 
#                   # pophaz = pm,
#                   surv.type = "surv.obs",
#                   # weights = list(agegr = c(0.2,0.4,0.4)),
#                   breaks = list(AGE = seq(0,100, 1)))
# st <- survtab(dt, print = NULL, #adjust = "agegr", 
#                   pophaz = pm,
#                   surv.type = "surv.rel",
#                   relsurv.method = "pp",
#                   # weights = list(agegr = c(0.2,0.4,0.4)),
#                   breaks = list(FUT = seq(0,5,1/12)))
# st <- survtab(dt, print = NULL, adjust = c("sex","agegr"), 
#                   pophaz = pm,
#                   surv.type = "surv.rel",
#                   relsurv.method = "pp",
#                   weights = list(sex = c(0.5, 0.5), agegr = c(0.2,0.4,0.4)),
#                   breaks = list(FUT = seq(0,5,1/12)))


# ag <- lexpand(sire, birth = "bi_date", entry = "dg_date", exit = "ex_date",
#               status = status %in% 1:2, pophaz = popmort, pp = TRUE,
#               fot = seq(0, 5, 1/12))
# pm2 <- copy(popEpi::popmort)
# setnames(pm2, c("year", "agegroup"), c("per", "age"))
# st <- survtab(ag, print = NULL, #adjust = c("sex","agegr"), 
#                   pophaz = pm2,
#                   surv.type = "surv.rel",
#                   relsurv.method = "pp",
#                   #weights = list(sex = c(0.5, 0.5), agegr = c(0.2,0.4,0.4)),
#                   breaks = list(fot = seq(0,5,1/12)))
detectEvents <- function(x, breaks, tol = .Machine$double.eps^0.5, by = "lex.id") {
  ## INTENTION: given a Lexis object, determines which rows
  ## have an event (a transition or end-point) within the window
  ## determined by breaks (a list of breaks as supplied to e.g. splitMulti).
  ## Usable with split and unsplit data, though it is best to do this
  ## before splitting for efficiency.
  ## NOTE: by should be a character vector specifying variables that identify
  ## unique subjects; the idea is that each subject only has one end-point
  ## and some transitions
  ## NOTE: breaks should be a list of breaks or NULL; if it is NULL,
  ## it is NOT checked whether observations were cut short by the breaks used.
  ## observations cut short are not any kind of events.
  ## OUTPUT: an integer vector coding events as follows:
  ## 0: no event within breaks (row cut short by breaks or subject 
  ##    has multiple rows, of which this is not an event)
  ## 1: transition within breaks
  ## 2: original end-point within breaks and no transition occured (i.e. censoring)
  if (!is.data.table(x)) stop("x must be a data.table; if you see this, send the package maintainer an email")
  # checkLexisData(x)
  if (!inherits(x, "Lexis")) stop("data not a Lexis object")
  if (!is.null(breaks)) {
    checkBreaksList(x, breaks)
    breaks[unlist(lapply(breaks, length)) == 0L] <- NULL
  }
  
  
  tmp <- list()
  oldKey <- key(x)
  if (length(oldKey) == 0L) {
    tmp$order <- makeTempVarName(x, pre = "order_")
    on.exit(setorderv(x, tmp$order), add = TRUE)
    on.exit(setcolsnull(x, tmp$order, soft = TRUE), add = TRUE)
    set(x, j = tmp$order, value = 1:nrow(x))
  } else on.exit(setkeyv(x, oldKey), add = TRUE)
    
  
  setkeyv(x, c(by, names(breaks)[1L]))
  setkeyv(x, by)
  ## rows that actually can be events: transitions and last rows by subject
  whTr <- x[, lex.Cst != lex.Xst]
  whLa <- !duplicated(x, fromLast = TRUE)
  whEv <- whTr | whLa
  
  if (!is.null(breaks)) {
    
    splitScales <- names(breaks)
    if (any(!splitScales %in% names(x))) stop("Following time scales missing from data that data was split by: ", paste0("'", setdiff(splitScales, names(x)), "'", collapse = ", "))
    
    brmax <- lapply(breaks, max)
    brmin <- lapply(breaks, min)
    
    ## detect rows residing within breaks window
    for (sc in splitScales) {
      z <- (x$lex.dur + x[[sc]])[whEv]
      tol_sc <- if (is.double(z)) tol else 0L
      
      ## NOTE: if max of orig values within breaks window, then all may be events
      if (!(max(z) + tol_sc < brmax[[sc]])) whEv[whEv] <- z < brmax[[sc]] - tol_sc
      if (!(min(z) - tol_sc > brmin[[sc]])) whEv[whEv] <- z > brmin[[sc]] + tol_sc
      
    }
    ## whEv now indicates rows that may be events AND which reside within breaks window. 
  }
  
  ## censorings are not transitions, but must reside within breaks window.
  whCe <- whLa & !whTr & whEv
  
  ## need to add event indicator to data since it has been reordered,
  ## reorder back old order, and return the event indicator.
  tmp$ind <- makeTempVarName(x, pre = "event_indicator_")
  on.exit(setcolsnull(x, delete = tmp$ind, soft = TRUE), add = TRUE)
  evInd <- as.integer(whEv)
  evInd <- ifelse(whCe, 2L, evInd)
  set(x, j = tmp$ind, value = evInd)
  
  if (length(oldKey) == 0L) {
    setkeyv(x, NULL)
    setorderv(x, tmp$order)
    set(x, j = tmp$order, value = NULL)
  } else setkeyv(x, oldKey)
  
  
  evInd <- x[[tmp$ind]]
  set(x, j = tmp$ind, value = NULL)
  on.exit(expr = {}, add = FALSE) ## removes on.exit expressions from earlier
  
  
  if (!identical(oldKey, key(x))) stop("keys do not match at function end; send an email to package maintainer if you see this")
  
  evInd
}









