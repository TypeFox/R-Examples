globalVariables(c("tab", "SE.A", "pp.table"))

comp.st.surv <- function(surv.var = "p", surv.expr = "1-d/(n.eff-n.cens/2)", 
                         SE.expr = "sqrt(p*(1-p)/(n.eff))", cumu = TRUE) {
  
  function(surv.table = tab, surv.by.vars = NULL) {
    
    tabname <- deparse(substitute(surv.table))
    #     tabname <- "surv.table"    
    by.expr <- NULL
    if (!is.null(surv.by.vars)) {
      surv.by.vars <- paste0("'", surv.by.vars, "'", collapse = " ,")
      by.expr <- paste0(", by = c( ", surv.by.vars, " )")
    }
    
    surv.expr <-  paste0(surv.var, " := ", surv.expr)
    tab.surv.expr <- paste0(tabname, "[, ", surv.expr,  by.expr, "]")
    #   surv.expr <- parse(text=surv.expr)
    #   surv.var <- parse(text=surv.var)
    SE.var <- paste0("SE.", surv.var)
    SE.expr <- paste0(SE.var, " := ", SE.expr)
    tab.SE.expr <- paste0(tabname, "[, ", SE.expr, ", ", by.expr, "]")
    #   SE.expr <- parse(text=SE.expr)
    #   SE.var <- parse(text=SE.var)
    cumuexpr <- paste0(tabname, "[, ", surv.var, " := cumprod(", surv.var, ")", by.expr, "]" )
    
    
    ## parent.frame(2): two steps upstream in environments
    pe <- function(obj, ..., env=parent.frame(2)) {
      eval(parse(text=obj), ..., envir=env)
    }
    
    pe(tab.surv.expr)
    
    
    if (cumu) pe(cumuexpr)
    
    ## standard error
    pe(tab.SE.expr)
    
    ## zero survival leads to zero SE
    minexpr <- paste0(tabname, "[",surv.var, "== 0, ", SE.var, " := 0]")
    pe(minexpr)
  }
}


comp.st.surv.obs.lif <- comp.st.surv(surv.var=  "surv.obs",
                                     surv.expr= "1-d/n.eff", 
                                     SE.expr=   "surv.obs*sqrt(  cumsum(  d/(n.eff*(n.eff-d))  )  )", ## old : sqrt(p*(1-p)/(n.eff))"
                                     cumu = TRUE)

comp.st.surv.obs.haz <- comp.st.surv(surv.var=  "surv.obs", 
                                     surv.expr= "exp(-delta*d/pyrs)", 
                                     SE.expr=   "surv.obs*sqrt(  cumsum(  delta^2*d/pyrs^2  )  )", ## old: sqrt(p*(1-p)/(n.eff))
                                     cumu = TRUE)


comp.st.r.e2.haz <- comp.st.surv(surv.var = "r.e2",
                                 surv.expr = "exp(-delta*(d-d.exp)/pyrs)",
                                 SE.expr = "SE.surv.obs/surv.exp",
                                 cumu=TRUE)

comp.st.r.e2.lif <- comp.st.surv(surv.var = "r.e2",
                                 surv.expr = "1-(d-d.exp)/(n.eff)",
                                 SE.expr = "SE.surv.obs/surv.exp",
                                 cumu=TRUE)

comp.st.r.pp.haz <- comp.st.surv(surv.var = "r.pp",
                                 surv.expr= "exp(-delta*(d.pp - d.exp.pp)/pyrs.pp)",
                                 SE.expr = "r.pp*sqrt(  cumsum(  delta^2*( d.pp.2)/ pyrs.pp^2  )  )",
                                 cumu=TRUE)

comp.st.r.pp.lif <- comp.st.surv(surv.var = "r.pp",
                                 surv.expr = "1-(d.pp-d.exp.pp)/(n.eff.pp)",
                                 SE.expr = "r.pp*sqrt(  cumsum(  d.pp.2/(n.eff^2)  )  )",
                                 cumu=TRUE)


## this function will calculate confidence intervals
## for obs & rel & net survivals
#' @import stats
comp.st.conf.ints <- function(tab = pp.table, al=0.05, surv="r.pp", transform ="log-log") {
  al <- 0.05
  zlo <- as.character(qnorm(al/2))
  zhi <- as.character(qnorm(1-al/2))
  SE.surv <- paste0("SE.",surv)
  surv.hi <- paste0(surv, ".hi")
  surv.lo <- paste0(surv, ".lo")
  
  pe <- function(...) {
    eval(parse(text=paste0(...)), envir=tab)
  }
  
  if (transform =="plain") {
    ## assume S(t)~N(mu, sigma)
    
    ex <- paste0(surv,  " ", zlo, "*", SE.surv)
    tab[, c(surv.lo)  := pe(ex)]
    ex <- paste0(surv, " +", zhi, "*", SE.surv)
    tab[, c(surv.hi)  := pe(ex)]
  }
  
  if (transform =="log-log") {
    ## assume log(H(t))~N(mu, sigma)
    ex <- paste0(SE.surv,"/(abs(log(",surv,"))*",surv,")")
    tab[,  SE.A := pe(ex)]
    
    ex <- paste0(surv, "^exp(", zhi, "*SE.A)")
    tab[, c(surv.lo)  := pe(ex)]
    ex <- paste0(surv, "^exp(", zlo, "*SE.A)")
    tab[, c(surv.hi)  := pe(ex)]
  }
  
  if (transform =="log") {
    ## assume log(S(t))~N(mu, sigma)
    ex <- paste0(SE.surv,"/",surv)
    tab[,  SE.A := pe(ex)]
    
    ex <- paste0(surv, "*exp(", zlo, "*SE.A)")
    tab[, c(surv.lo)  := pe(ex)]
    ex <- paste0(surv, "*exp(", zhi, "*SE.A)")
    tab[, c(surv.hi)  := pe(ex)]
  }
  
  ## zero SE means zero uncertainty means lo=hi=estimate
  tab[get(SE.surv) == 0, c(surv.lo, surv.hi) := get(surv)]
  
  tab[]
}



# x <- Lexis(data=sire[1,], entry = list(fot=0, per=get.yrs(dg_date), age=dg_age), 
#            exit=list(per=get.yrs(ex_date)), exit.status=status)
# x <- splitMulti(x, breaks = list(fot=seq(0, 5, by = 1/12), per=1994:2013, age = 0:150))
# x[, surv.int := cut(fot, seq(0, 5, 1/12) - .Machine$double.eps^0.5, labels = FALSE)]
# x <- cutLowMerge(x, popmort, by.x = c("sex","per", "age"), 
#                  by.y = c("sex", "year", "agegroup"), 
#                  mid.scales = c("per", "age"), all.x = TRUE, all.y = FALSE)
# comp_pp_weights(x, surv.scale = "fot", breaks = seq(0, 5, 1/12), haz = "haz", style = "delta")

comp_pp_weights <- function(lex, surv.scale = "fot", breaks = NULL, haz = "haz", style = "delta", verbose = FALSE) {
  ppTime <- proc.time()
  ## input: a split Lexis object (data.table) and the time scale to compute
  ## pp weights over; lex must also contain 'haz', the population
  ## (expected) hazard level for each row
  TF <- environment()
  
  lex.id <- pp <- NULL ## APPEASE R CMD CHECK
  
  style <- match.arg(style, c("delta", "actual"))
  if (!is.data.table(lex)) stop("lex must be a data.table")
  
  all_names_present(lex, c(haz, surv.scale, "lex.id", "lex.dur"))
  
  if ("pp" %in% names(lex)) stop("Variable named 'pp' existed in data when attempting to compute pohar-perme weights; 'pp' is reserved so you should delete or rename that variable")
  
  
  if (!identical(key(lex), c("lex.id", surv.scale))) stop("lex must be a data.table keyed by lex.id and the survival time scale")
  
  
  breaks <- sort(unique(breaks)) - .Machine$double.eps^0.5
  
  ## need a bunch of temporary variable names to compute pp weights
  ## inside the data set without overwriting anything existing.
  ## this will take about 0.1 secs.
  tmpSI <- makeTempVarName(data = lex, pre = "surv.int_")
  tmpSIstart <- makeTempVarName(data = lex, pre = "surv.int.start_")
  tmpSIstop <- makeTempVarName(data = lex, pre = "surv.int.stop_")
  tmpSIlength <- makeTempVarName(data = lex, pre = "surv.int.length_")
  on.exit(setcolsnull(lex, delete = c(tmpSI, tmpSIstart, tmpSIstop, tmpSIlength, tmpPS, tmpPCS, tmpPCSM)), add = TRUE)
  
  set(lex, j = tmpSI, value = cut(lex[[surv.scale]], breaks, labels = FALSE))
  set(lex, j = tmpSIstop, value = breaks[-1][lex[[tmpSI]]])
  set(lex, j = tmpSIstart, value = breaks[-length(breaks)][lex[[tmpSI]]])
  set(lex, j = tmpSIlength, value = lex[[tmpSIstop]] - lex[[tmpSIstart]])
  
  tmpPS <- makeTempVarName(data = lex, pre = "pop.surv_")
  tmpPCS <- makeTempVarName(data = lex, pre = "pop.cumsurv_")
  tmpPCSM <- makeTempVarName(data = lex, pre = "pop.cumsurv.mid_")

  ## conditional survs
  if (verbose) condSurvTime <- proc.time()
  
  set(lex, j = tmpPS, value = exp(-lex[[haz]]*lex$lex.dur))
  ## till end of each interval...
  lex[, c(tmpPCS) := list(exp(-cumsum(.SD[[1L]]*lex.dur))),
      by = lex.id, .SDcols = c(haz, "lex.dur")]
  ## till start of each interval
  set(lex, j = tmpPCS, value = lex[[tmpPCS]] / lex[[tmpPS]]) 
  if (verbose) cat("Time taken by computing expected survivals up to start of each interval for each lex.id: ", timetaken(condSurvTime), "\n")
  
  ## pohar-perme weighting by expected cumulative survival. approximation:
  ## cumulative survival up to either middle of remaining surv.int (not individual-specific)
  ## or up to middle of subject's follow-up in each row (individual-specific)
  ## difference: e.g. 2 rows within a surv.int have either the same or different pp-weights
  if (style == "actual") {
    set(lex, j = tmpPCSM, value = lex[[tmpPCS]] * (lex[[tmpPS]])^0.5)
  }
  if (style == "delta") {
    if (verbose) deltaTime <- proc.time()
    
    setkeyv(lex, c("lex.id", tmpSI))
    ## expected survival up to middle of remaining time in surv.int
    ## cumulation starting from first record for subject in each surv.int
    
    ## some records are the only one for a lex.id in a surv.int; these are easy
    first_in_surv.int <- !duplicated(lex, fromLast = FALSE)
    last_in_surv.int <- !duplicated(lex, fromLast = TRUE)
    only_in_surv.int <- first_in_surv.int & last_in_surv.int
    
    lex[only_in_surv.int, c(tmpPCSM) := .SD[[1L]] * exp(.SD[[2L]] * (.SD[[3L]] - .SD[[4L]])/2L),
        .SDcols = c(tmpPCS, haz, tmpSIstop, surv.scale)]
    
    ## more complicated with many records in a surv.int per lex.id
    if (any(!only_in_surv.int)) {
      
      tmpSImid <- makeTempVarName(lex, pre = "surv.int.mid_")
      dist <- makeTempVarName(lex, pre = "dist_")
      on.exit(setcolsnull(lex, delete = c(dist, tmpSImid)), add = TRUE)
      
      ## middle point of survival interval
      set(lex, j = tmpSImid, value = (lex[[tmpSIstop]] - lex[[tmpSIstart]])/2L)
      
      ## distance from remaining surv.int mid-point starting from start of 
      ## record; or at most lex.dur; for integration
      set(lex, j = dist, value = pmin(lex$lex.dur, lex[[tmpSImid]]))
      ## some records after mid-point can have negative fot.dist at this point
      set(lex, j = dist, value = pmax(lex[[dist]], 0))
      
      ## some lex.id are censored / die before mid of surv.int; last record
      ## must reach its fot.dist at least up to the mid (or be zero due to above)
      lex[last_in_surv.int, c(dist) := pmax((.SD[[1L]] - .SD[[2L]])/2L, 0), 
          .SDcols = c(tmpSIstop, surv.scale)]
      
      byTime <- proc.time()
      
      ## from start of first in surv.int till mid point
      ## step by step for efficiency...
      lex[!only_in_surv.int, c(tmpPCSM) := .SD[[1L]] * .SD[[2L]], .SDcols = c(haz, dist)]
      lex[!only_in_surv.int, c(tmpPCSM) := lapply(.SD, sum), 
          .SDcols = c(tmpPCSM), by = c("lex.id", tmpSI)]
      lex[!only_in_surv.int, c(tmpPCS) := .SD[1L], 
          by = c("lex.id", tmpSI), .SDcols = c(tmpPCS)]
      lex[!only_in_surv.int, c(tmpPCSM) := .SD[[1L]] * exp(-.SD[[2L]]), 
          .SDcols = c(tmpPCS, tmpPCSM)]
      
      ## todo: alternate faster method for integration!
      setcolsnull(lex, delete = c(dist))
      if (verbose) cat("Time taken by extra computation due to style 'delta': ", timetaken(deltaTime), "\n")
    }
    
    rm(first_in_surv.int, last_in_surv.int, only_in_surv.int)
    if (verbose) cat("Time taken by computation of Pohar-Perme weights: ", timetaken(ppTime), "\n")
  }
  
  setkeyv(lex, c("lex.id", surv.scale))
  
  lex[, pp := 1/.SD, .SDcols = tmpPCSM]
  
  invisible(lex[])
}

comp_pp_weighted_figures <- function(lex, haz = "haz", pp = "pp", event.ind = NULL, by = "lex.id") {
  ## PURPOSE: given a split Lexis object with a column of pre-computed
  ## pohar-perme weights, computes pp-weighted:
  ## * person-time (lex.dur*pp)
  ## * at-risk indicator
  ## * event counts
  ## * event counts multiplied with pp^2
  ## * expected event counts
  ## events are transitions and end points as detected by detectEvents,
  ## and include censorings.
  ## OUTPUT: a DT of pp-weighted things.
  
  checkLexisData(lex, check.breaks = TRUE)
  if (!is.data.table(lex)) stop("lex must be a data.table")
  
  all_names_present(lex, c(pp, haz, event.ind, by))
  
  ## NOTE: we want to avoid conflicts with possible variable names in lex
  ## (e.g. an existing column named ptime.pp might conflict with ptime.pp char vec)
  e <- environment()
  
  if (is.null(event.ind)) {
    event.ind <- makeTempVarName(lex, pre = "event_indicator_")
    on.exit(setcolsnull(lex, delete = event.ind, soft = TRUE), add = TRUE)
    set(lex, j = event.ind, value = detectEvents(lex, breaks = attr(lex, "breaks"), by = by))
    all_names_present(lex, event.ind)
  }
  ## data.table is probably faster in this than simply using vectors
  idt <- data.table(1:2)
  names(idt) <- event.ind
  haveEvents <- sort(lex[idt, on = event.ind, which = TRUE])
  evtab <- lex[haveEvents, .(obs = .N), by = list(lex.Cst, lex.Xst)]
  set(evtab, j = "obs", value = NULL)
  
  events <- paste0("from", evtab$lex.Cst, "to", evtab$lex.Xst)
  ppVars <- c(paste0(events, ".pp"), "ptime.pp", "d.exp.pp")
  
  
#   ## build table to join with lex to limit to rows with events -----------------
#   ## non-event-rows have zero valued pp-weighted event counts, naturally
#   set(evtab, j = event.ind, value = NULL)
#   evtab <- rbindlist(list(evtab, evtab))
#   set(evtab, j = event.ind, value = rep(1:2, each = nrow(evtab)/2L))
#   setkeyv(evtab, c("lex.Cst", "lex.Xst", event.ind))
  
  ## pp-weighted events --------------------------------------------------------
  ## NOTE: a pp-weighted event is simply the pp weight where the event occurs
  ## and zero otherwise (0L/1L times pp-weight)
  
  evN <- length(events)
  evdt <- data.table(rn = rep(1:nrow(lex), times = evN))
  set(evdt, j = "eventType", value = factor(rep(1:evN, each = nrow(lex)), levels = 1:evN, labels = events))
  set(evdt, j = "pp", value = lex[[pp]])
  
  
  ## need to still determine which rows are not their eventType's events -------
  
  rnVar <- makeTempVarName(lex, pre = "rowNum_")
  on.exit(setcolsnull(lex, delete = rnVar, soft = TRUE), add = TRUE)
  set(lex, j = rnVar, value = 1:nrow(lex))
  
  ## row numbers in lex by event type
  rowNums <- lex[haveEvents][evtab, .(rowNum = .SD[[e$rnVar]]), by = .EACHI, on = c("lex.Cst", "lex.Xst"), .SDcols = rnVar]
  ## multiply to accommodate expanded evdt data
  
  rowNum <- NULL ## appease R CMD CHECK
  rowNums[, rowNum := rowNum+(.GRP-1L)*nrow(lex), by = list(lex.Cst, lex.Xst)]
  noEvents <- setdiff(1:nrow(evdt), rowNums$rowNum)
  
  set(evdt, i = noEvents, j = "pp", value = 0)
  
  evdt <- dcast.data.table(evdt, rn ~ eventType, value.var = "pp")
  
  setorderv(evdt, "rn")
  set(evdt, j = "rn", value = NULL)
  
  ## ptime.pp & d.exp.pp -------------------------------------------------------
  set(evdt, j = "ptime.pp", value = lex$lex.dur * lex[[pp]])
  set(evdt, j = "d.exp.pp", value = lex$lex.dur * lex[[haz]] * lex[[pp]])
  
  setcolorder(evdt, c("ptime.pp", "d.exp.pp", sort(events)))
  setnames(evdt, events, paste0(events, ".pp"))
  
  ## pp-weighted at-risk indicators --------------------------------------------
  ## these will be n.pp at the aggregate level.
  set(evdt, j = "at.risk.pp", value = lex[[pp]]*1L)
  
  ## events multiplied with pp-weight again ------------------------------------
  ## (i.e. event times pp squared)
  
  set(evdt, j = paste0(events, ".pp.2"), 
      value = evdt[, .SD, .SDcols = paste0(events, ".pp")]*lex[[pp]])
  
  return(evdt[])
  
}



test_empty_surv_ints <- function(x, by = NULL, sum.over = NULL, test.var = "pyrs", show.by = NULL) {
  
  x <- copy(x)
  oc <- class(x)
  setDT(x)
  setattr(x, "class", oc)
  
  surv.int <- NULL ## APPEASE R CMD CHECK
  
  all_names_present(x, by, msg = "Missing variable(s) %%VARS%% from data when inspected for empty survival intervals. If you see this, send a message to the package maintainer.")
  all_names_present(x, sum.over, msg = "Missing variable(s) %%VARS%% from data when inspected for empty survival intervals. If you see this, send a message to the package maintainer.")
  all_names_present(x, test.var, msg = "Missing variable(s) %%VARS%% from data when inspected for empty survival intervals. If you see this, send a message to the package maintainer.")
  
  if (any(!sum.over %in% by)) stop("sum.over must be a subset of by.")
  if (length(show.by) == 0L) show.by <- by
  if (length(by) != length(show.by)) {
    stop("Internal error: length(sum.over) != length(show.sum.over). ",
         "If you see this, complain to the package maintainer.")
  }
  
  wh_sum.to <- !by %in% sum.over
  sum.to <- by[wh_sum.to]
  show.sum.to <- show.by[wh_sum.to]
  
  if (length(sum.to) == 0L) sum.to <- show.sum.to <- NULL
  
  tmpTV <- makeTempVarName(x, pre = "testValues_")
  tmpDiff <- makeTempVarName(x, pre = "diff_")
  on.exit(setcolsnull(x, delete = c(tmpTV, tmpDiff)))
  
  ## first check empty surv.ints are all consecutive...
  ## consecutiveness: e.g. out of 10 surv ints, 6-10 are empty.
  ## non-consecutiveness: e.g. out of 10 surv ints, 6-8 are empty.
  ## (then 9-10 will have NA estimates as well.)
  setkeyv(x, c(by, "surv.int"))
  
  ct <- x[, lapply(.SD, sum), keyby = c(sum.to, "surv.int"), .SDcols = test.var]
  setnames(ct, length(ct), tmpTV)
  ct <- ct[ct[[tmpTV]] > 0L, diff(surv.int), keyby = eval(sum.to)]
  setnames(ct, length(ct), tmpDiff)
  ct <- ct[ct[[tmpDiff]] > 1L]
  ## THE IDEA: if the difference in the number of the survival interval
  ## is > 1, it means there is at least one row between two non-empty
  ## intervals, i.e. non-consecutively.
  
  ## we keep non-consecutively bad surv.int stratas in entirety for inspection
  if (nrow(ct) > 0L) {
    msg <- paste0("The total person-time was zero in some survival intervals")
    
    if (!is.null(sum.to)) {
      msg <- paste0(msg, ", when summed to the variable(s) ", 
                    paste0("'", show.sum.to, "'", collapse = ", "),
                    " (i.e. over all other variables, if any)") 
    } else {
      msg <- paste0(msg, " summed to the margins (over any stratifying ",
                    "/ adjusting variables)")
    }
      
    msg <- paste0(msg, " _non-consecutively_, i.e. some intervals after an ",
                  "empty interval had person-time in them. ",
                  "Keeping all survival ",
                  "intervals with some estimates as NA for inspection.")
    message(msg)
  } else {
    ## we leave out intervals that are empty consecutively (e.g. all from 5-10)
    x[, c(tmpTV) := lapply(.SD, sum), by=c(sum.to, "surv.int"), .SDcols = test.var]
    x <- x[x[[tmpTV]] > 0L]
    setcolsnull(x, tmpTV)
  }
  
  x
}




comp_e1 <- function(x, breaks, pophaz, survScale, by = NULL, id = "lex.id", immortal = TRUE, verbose = FALSE) {
  ## INTENTION: given a Lexis data set x,
  ## computes Ederer I expected survival curves till end of follow-up
  ## by 'by' unless individual = TRUE.
  TF <- environment()
  
  haz <- surv.exp <- NULL # R CMD CHECK appeasement
  ## check ---------------------------------------------------------------------
  checkLexisData(x)
  checkBreaksList(x, breaks)
  ph <- data.table(pophaz)
  tmpHaz <- makeTempVarName(x, pre = "pop_haz_")
  setnames(ph, "haz", tmpHaz)
  checkPophaz(x, ph, haz.name = tmpHaz)
  
  byErr <- paste0("Internal error (probably): work data did not have ",
                  "variable(s) %%VARS%%. If your supplied data has them, ",
                  "complain to the package maintainer.")
  all_names_present(x, c(by, id, survScale), msg = byErr)
  
  if (length(id) != 1L) {
    stop("Argument id must be of length 1.")
  }
  ## split ---------------------------------------------------------------------
  pt <- proc.time()
  oldBreaks <- attr(x, "breaks")
  allScales <- attr(x, "time.scales")
  if (!survScale %in% allScales) {
    stop("survScale '", survScale, "' not a time scale in the Lexis object. ",
         "(Possibly internal error - ensure you have that time scale in data. ",
         "If not, complain to the package maintainer.")
  }
  keepVars <- c(allScales, "lex.dur", "lex.id", "lex.Cst", 
                "lex.Xst", by, id, setdiff(names(ph), tmpHaz))
  keepVars <- unique(keepVars)
  y <- subsetDTorDF(x, select = keepVars)
  y <- setDT(copy(y))
  forceLexisDT(y, breaks = oldBreaks, allScales = allScales, key = TRUE)
  
  ## won't use statuses for anything
  y[, c("lex.Cst", "lex.Xst") := NULL]
  y[, c("lex.Cst", "lex.Xst") := 0L]
  
  if (immortal) {
    ## set lex.dur to infinite. this assumes that the subjects never leave
    ## follow-up (which Ederer I assumes)
    storage.mode(y$lex.dur) <- "double"
    id_last <- !duplicated(y, by = "lex.id", fromLast = TRUE)
    y[TF$id_last, lex.dur := Inf]
  }
  
  y <- intelliCrop(y, breaks = breaks, allScales = allScales)
  y <- intelliDrop(y, breaks = breaks)
  setDT(y)
  forceLexisDT(y, breaks = oldBreaks, allScales = allScales, key = TRUE)
  setkeyv(y, c(id, survScale))
  
  y <- splitMulti(y, breaks = breaks, drop = FALSE, merge = TRUE)
  
  if (verbose) cat("Time taken by splitting: ", timetaken(pt), ".\n", sep = "")
  
  ## merge pop haz -------------------------------------------------------------
  pt <- proc.time()
  mergeVars <- intersect(names(y), names(ph))
  mergeScales <- intersect(allScales, mergeVars)
  if (length(mergeScales) == 0L) mergeScales <- NULL
  mergeCats <- setdiff(mergeVars, mergeScales)
  if (length(mergeCats) == 0L) mergeCats <- NULL
  
  y <- cutLowMerge(y, ph, by = mergeVars, 
                   mid.scales = mergeScales, old.nums = TRUE,
                   all.x = TRUE, all.y = FALSE)
  setDT(y)
  if (verbose) cat("Time taken by merging pophaz: ", timetaken(pt), ".\n", sep = "")
  
  ## ederer I computation ------------------------------------------------------
  ## prep temp surv.int var
  tmpSI <- makeTempVarName(x, pre = "surv_int_")
  setkeyv(y, c(by, id, allScales[1L]))
  y[, c(tmpSI) := cut(y[[survScale]], breaks[[survScale]],
                      right=FALSE,labels=FALSE)]
  setkeyv(y, c(by, tmpSI))
  
  ## EDERER I: integral of the weighted average expected hazard,
  ## where the weights are the subject-specific expected survival
  ## probabilities.
  ## 1) compute integral of hazard over an interval t_i by id
  ##    (NOT cumulative hazard from zero till end of interval t_i)
  ##    This sums over multiple rows a subject may have within one
  ##    survival interval due to splitting by multiple time scales.
  pt <- proc.time()
  set(y, j = tmpHaz, value = y[[tmpHaz]]*y$lex.dur)
  y <- y[, lapply(.SD, sum), keyby = eval(unique(c(by, id, tmpSI))), 
         .SDcols = c(tmpHaz)]
  setnames(y, ncol(y), tmpHaz)
  
  ## reverse temp names - need to be able to refer to haz without temp var
  ## to enable correct cumsum() below
  avoid <- c(names(y), tmpHaz, "surv.exp")
  tmpBy <- makeTempVarName(names = avoid, pre = by)
  tmpID <- makeTempVarName(names = c(avoid, tmpBy), pre = id)
  if (length(by)) {
    setnames(y, by, tmpBy)
    if (id %in% by) tmpID <- id <- tmpBy[by == id]
  } else {
    tmpBy <- by <- NULL
  }
  setnames(y, id, tmpID)
  setnames(y, tmpHaz, "haz")
  if (verbose) cat("Time taken by 1): ", timetaken(pt), ".\n", sep = "")
  
  ## 2) expected cum.haz. over intervals t_1 -> t_i by id...
  ##   (no cumulative exp.surv yet)
  pt <- proc.time()
  y[, surv.exp := cumsum(haz), by = eval(tmpID)]
  
  if (verbose) cat("Time taken by 2): ", timetaken(pt), ".\n", sep = "")
  ## 3) cumulative surv.exp till end of interval t_i by id...
  pt <- proc.time()
  y[, surv.exp := exp(-surv.exp)]
  
  if (verbose) cat("Time taken by 3): ", timetaken(pt), ".\n", sep = "")
  
  ## 4) The Ederer I expected (marginal) survivals for intervals t_i 
  pt <- proc.time()
  y <- y[, .(surv.exp = mean(surv.exp)), by = eval(c(tmpBy, tmpSI))]
  if (verbose) cat("Time taken by 4): ", timetaken(pt), ".\n", sep = "")
  
  if ("surv.exp" %in% by) {
    by[by == "surv.exp"] <- makeTempVarName(y, pre = "surv.exp")
  }
  if (length(by)) setnames(y, tmpBy, by)
  
  setcolorder(y, c(by, tmpSI, "surv.exp"))
  setkeyv(y, c(by, tmpSI))
  setnames(y, tmpSI, survScale)
  br <- breaks[[survScale]]
  br <- br[-1]
  y[, c(survScale) := br[y[[survScale]]]]
  
  
  y[]
}




detectSurvivalTimeScale <- function(lex, values) {
  
  checkLexisData(lex)
  
  allScales <- attr(lex, "time.scales")
  
  allScVals <- lapply(allScales, function(ch) lex[[ch]])
  names(allScVals) <- allScales
  
  whSurvScale <- lapply(allScVals, function(col) {
    identical(col, values)
  })
  whSurvScale <- unlist(whSurvScale)
  
  if (sum(whSurvScale) == 0L) {
    whSurvScale <- lapply(allScVals, function(col) {
      isTRUE({
        all.equal(col, values, scale = 1L, 
                  check.attributes = FALSE,
                  tolerance = .Machine$double.eps ^ 0.5)
      })
    })
    whSurvScale <- unlist(whSurvScale)
  }
  if (sum(whSurvScale) == 0L) {
    
    dt <- as.data.table(allScVals)
    dt <- cbind(dt, data.table(testValues = values))
    on.exit(print(dt))
    
    stop("Could not determine which time scale was used. The formula MUST ",
         "include the time scale used within a Surv() call (or a Surv object),",
         " e.g. Surv(FUT, lex.Xst) ~ sex. Note that the 'time' argument is ",
         "effectively (and exceptionally) used here to denote the times at ",
         "the beginning of follow-up to identify the time scale existing in ",
         "the supplied data to use. If you are sure you are mentioning a ",
         "time scale in the formula in this manner, complain to the ",
         "package maintainer. The table printed below contains the time ",
         "scales tested against and the values that were supplied as the last ",
         "column.")
  }
  survScale <- allScales[whSurvScale]
  survScale
  
}
