#############################################################################
##
##  Stage 3: Approximation Errors for Numerical Inversion
##
#############################################################################

## --------------------------------------------------------------------------
## Create page for computing approximation errors of numerical inversion

experim.inverror <- function (main, win, nb) {

  ## get data about generator object
  method.str <- tag(main,"method")$symbol

  ## currently it is only implemented for method PINV
  if(method.str!="PINV") return()
      
  ## add page for tests
  paerror <- ggroup(horizontal=FALSE, spacing=10, container=nb, label="Approximation Error")

  ## mark page as protected (i.e., cannot be closed)
  tag(win,"protected") <- append(tag(win,"protected"), length(nb))

  ## add functions for toolbar handlers.
  ## functions must be of type 'function(win)'
  pageno <- as.character(length(nb))
  tag(win,"run")[[pageno]] <- experim.ierror
  tag(win,"plot")[[pageno]] <- experim.ierror.plot

  ## Distribution
  params.distr.gwl <- experim.aerror.distr(main, win, paerror)

  ## we group frames for histogram and tests horizontally
  aerror.grp <- ggroup(horizontal=TRUE, expand=TRUE, spacing=15, container=paerror)
  
  ## Data (type of error)Histogram
  params.aerror.gwl <- experim.aerror.data(main, win, aerror.grp)

  ## Show result
  params.result.gwl <- experim.aerror.result(main, win, aerror.grp)

  ## store list of widgets in window
  tag(win,"params.aerror.gwl") <- c(params.distr.gwl, params.aerror.gwl, params.result.gwl)

  ## return page
  return (paerror)
}


## --------------------------------------------------------------------------
## Approximation error: Distribution 

experim.aerror.distr <- function (main, win, group) {

  ## similar to experim.hist.distr
  params.gwl <- experim.hist.distr(main, win, group)

  svalue(params.gwl[["cdf.cbx"]])  <- TRUE
  svalue(params.gwl[["icdf.cbx"]]) <- FALSE

  enabled(params.gwl[["cdf.cbx"]])  <- FALSE
  enabled(params.gwl[["icdf.cbx"]]) <- FALSE
         
  ## return list of widgets
  return (params.gwl)
}


## --------------------------------------------------------------------------
## Approximation error: Data

experim.aerror.data <- function (main, win, group) {
  
  ## store widgets for parameters in a list
  params.gwl <- list()

  ## create a frame for inserting parameters for histogram
  frame <- gframe(" Type of Error ", spacing=10, horizontal=FALSE, container=group)
  addSpace(frame,10)
  
  ## type of error
  type.rbx <- gradio(c("u-error   ","absolute x-error   ","relative x-error"),
                     horizontal=TRUE, container=frame,
                     handler=function(h,...){experim.aerror.update(win)})
  params.gwl[["type.rbx"]] <- type.rbx
  
  glabel("Remark: x-errors are not well suited for \nmeasuring errors in numerical inversion",
         container=frame)

  ## create a table for the layout of widgets
  tbl <- glayout(container=frame)

  ## sample size
  n.sample.edt   <- gedit("10000", coerce.with=as.integer, width=10, container=tbl)
  tbl[1,2:3,anchor=c(-1,0)] <- "sample size"
  tbl[1,4,anchor=c(-1,0)] <- params.gwl[["n.sample.edt"]] <- n.sample.edt
  tooltip(n.sample.edt) <- "Size of random sample"

  ## maximal tolerated error
  ures.edt <- tag(main,"method.params.gwl")[["uresolution"]]
  tol <- if(is.null(ures.edt)) 1.e-10 else svalue(ures.edt)
  tol.edt <- gedit(tol, coerce.with=as.numeric, width=10, container=tbl)
  tbl[2,2:3,anchor=c(-1,0)] <- "maximal tolerated error"
  tbl[2,4,anchor=c(-1,0)] <- params.gwl[["tol.edt"]] <- tol.edt
  tooltip(tol.edt) <- "Maximal tolerated error. Also used to scale the plot. Optional"
  
  ## u-domain
  lb.edt   <- gedit("0", coerce.with=as.numeric, width=10, container=tbl)
  ub.edt   <- gedit("1", coerce.with=as.numeric, width=10, container=tbl)
  tbl[3,2,anchor=c(-1,0)] <- "domain in u-scale"
  tbl[3,3,anchor=c(+1,0)] <- "["
  tbl[3,4,anchor=c(-1,0)] <- params.gwl[["lb.edt"]] <- lb.edt
  tbl[3,5,anchor=c(-1,0)] <- ","
  tbl[3,6,anchor=c(-1,0)] <- params.gwl[["ub.edt"]] <- ub.edt
  tbl[3,7,anchor=c(-1,0)] <- "]"
  tooltip(lb.edt) <- "Lower boundary for domain of investigation for inverse CDF"
  tooltip(ub.edt) <- "Upper boundary for domain of investigation for inverse CDF"
  
  ## return list of widgets
  return (params.gwl)
}


## --------------------------------------------------------------------------
## Approximation error: Result

experim.aerror.result <- function(main, win, group) {

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## create a frame for displaying results of GoF tests
  frame <- gframe(" Approximaton Errors ", spacing=10, horizontal=FALSE, container=group)
  
  ## create a table for the layout of widgets
  tbl <- glayout(container=frame)

  ## print run time
  tbl[2,2,anchor=c(-1,0)] <- "Total run time:"
  tbl[2,3,anchor=c(-1,0)] <- params.gwl[["runtime.lbl"]] <- glabel("- - : - - . - - -",container=tbl)
  tbl[2,4] <- "  "
  
  ## print simple statistics on error
  tbl[4,2,anchor=c(-1,0)] <- "maximum error:"
  tbl[4,3,anchor=c(-1,0)] <- params.gwl[["aerror.max.lbl"]] <- glabel("--",container=tbl)
  
  tbl[5,2,anchor=c(-1,0)] <- "mean absolute error:"
  tbl[5,3,anchor=c(-1,0)] <- params.gwl[["aerror.mad.lbl"]] <- glabel("--",container=tbl)
  
  tbl[6,2,anchor=c(-1,0)] <- "root of MSE:"
  tbl[6,3,anchor=c(-1,0)] <- params.gwl[["aerror.rmse.lbl"]] <- glabel("--",container=tbl)
  
  ## return list of widgets
  return (params.gwl)
}


## --------------------------------------------------------------------------
## update handler

experim.aerror.update <- function (win) {

  ## get data
  params.gwl <- tag(win,"params.aerror.gwl")

  if (svalue(params.gwl[["type.rbx"]], index=TRUE) == 1) {
    ## u-error
    svalue(params.gwl[["cdf.cbx"]])  <- TRUE
    svalue(params.gwl[["icdf.cbx"]]) <- FALSE
  }
  else {
    ## x-error
    svalue(params.gwl[["cdf.cbx"]])  <- FALSE
    svalue(params.gwl[["icdf.cbx"]]) <- TRUE
  }
}


## --------------------------------------------------------------------------
## Compute RVGT frequency table and compute p-values for tests

experim.ierror <- function(win) {

  ## update statusbar
  ## sb <- tag(win,"sb")
  ## svalue(sb) <- "Computing histogram. This may take some time ..."

  ## list of widgets with data
  params.gwl <- tag(win,"params.aerror.gwl")

  ## sample size
  n <- svalue(params.gwl[["n.sample.edt"]])
  
  if (n<100) {
    error.message("sample size < 100")
    tag(win,"ierror") <- NULL
    return(NULL)
  }

  ## domain for inverse CDF
  udom <- c(svalue(params.gwl[["lb.edt"]]), svalue(params.gwl[["ub.edt"]]))
  if (! (udom[1] >= 0. && udom[2] <= 1. && udom[1] < udom[2])) {
    error.message("domain is missing or violates 0 <= lb < ub <= 1")
    tag(win,"ierror") <- NULL
    return(NULL)
  }

  ## CDF and inverse CDF
  cdf.str  <- if( isTRUE(svalue(params.gwl[["cdf.cbx"]])))  svalue(params.gwl[["cdf.edt"]])  else NULL
  icdf.str <- if( isTRUE(svalue(params.gwl[["icdf.cbx"]]))) svalue(params.gwl[["icdf.edt"]]) else NULL

  if (is.null(cdf.str) && is.null(icdf.str)) {
    error.message("CDF / quantile function required")
    tag(win,"ierror") <- NULL
    return(NULL)
  }
  
  pdist <- tryCatch(eval(parse(text=paste(cdf.str))),
                    error=function(e){
                      error.message(e$message, title="UNU.RAN - Error: CDF"); NULL })

  qdist <- tryCatch(eval(parse(text=paste(icdf.str))),
                    error=function(e){
                      error.message(e$message, title="UNU.RAN - Error: quantile"); NULL })

  if (((svalue(params.gwl[["type.rbx"]], index=TRUE) == 1) && is.null(pdist)) ||
      ((svalue(params.gwl[["type.rbx"]], index=TRUE) != 1) && is.null(qdist)) ) {
    tag(win,"ierror") <- NULL
    return(NULL)
  }

  
  ## truncated distribution
  if (isTRUE(svalue(params.gwl[["istruncated"]]))) {
    trunc <- params.gwl[["trunc.val"]]
  }
  else {
    trunc <- NULL
  }

  ## numerical inversion method
  aqdist <- function(u) { uq(tag(win,"gen"), u) }

  ## start timer
  time.start <- proc.time()
  
  ## estimate approximation errors
  ierror <-
    tryCatch(switch(svalue(params.gwl[["type.rbx"]], index=TRUE),
                    "1" = { uerror(n, aqdist, pdist, trunc=trunc, udomain=udom) },
                    "2" = { xerror(n, aqdist, qdist, trunc=trunc, udomain=udom, kind="abs") },
                    "3" = { xerror(n, aqdist, qdist, trunc=trunc, udomain=udom, kind="rel") }),
             error=function(e){
               error.message(e$message, title="UNU.RAN - Error: Cannot compute approximation error")
               NULL
             } )

  if (is.null(ierror)) {
    tag(win,"ierror") <- NULL
    return(NULL)
  }

  ## stop timer
  run.time <- (proc.time() - time.start)[3]  ## "elapsed" time

  ## print result into frame
  svalue(params.gwl[["aerror.max.lbl"]])  <- signif(max(ierror$max),3)
  svalue(params.gwl[["aerror.mad.lbl"]])  <- signif(sum(ierror$mad)/ierror$res,3)
  svalue(params.gwl[["aerror.rmse.lbl"]]) <- signif(sqrt(sum(ierror$mse)/ierror$res),3)

  ## print run time into frame
  min <- as.integer(floor(run.time / 60))
  sec <- as.integer(floor(run.time - min*60))
  msec <- as.integer(round(1000*(run.time - min*60 - sec)))
  svalue(params.gwl[["runtime.lbl"]]) <- sprintf("%02d:%02d.%03d", min, sec, msec)
  
  ## store table in window and reset counter for plots
  tag(win,"ierror") <- ierror
  tag(win,"ierror.plotted") <- FALSE

  ## update statusbar
  ## svalue(sb) <- "Approximation errors estimated."

  ## return table
  return (ierror)
}                          
  

## --------------------------------------------------------------------------
## plot approximation errors

experim.ierror.plot <- function(win) {

  ## maximal accepted error
  tol <- svalue(tag(win,"params.aerror.gwl")[["tol.edt"]])
  if (isTRUE(tol<=0.)) {
    error.message("maximal tolerated error must be empty or > 0")
    return()
  }

  ## get table of errors
  ierror <- tag(win,"ierror")
  if (is.null(ierror) || isTRUE(tag(win,"ierror.plotted"))) {
    ## recompute
    ierror <- experim.ierror(win)
  }
  if (is.null(ierror)) return()

  ## page label
  label <- switch(svalue(tag(win,"params.aerror.gwl")[["type.rbx"]], index=TRUE),
                  "1" = "u-error",
                  "2" = "abs x-error",
                  "3" = "rel x-error" )

  ## approximation errors
  ggraphics(container=tag(win,"nb"), label=label)
  Sys.sleep(0.5)
  plot(ierror,tol=tol)
  tag(win,"ierror.plotted") <- TRUE
}

## --------------------------------------------------------------------------
