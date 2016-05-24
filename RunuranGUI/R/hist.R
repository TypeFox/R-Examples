#############################################################################
##
##  Stage3: Histogram Based Tests
##
#############################################################################

## --------------------------------------------------------------------------
## Create page for histogram based tests

experim.histtest <- function (main, win, nb) {

  ## check type of distribution
  if (tag(main,"distribution")$type != "continuous")
    return(NULL)
  
  ## add page for tests
  ptests <- ggroup(horizontal=FALSE, spacing=10, container=nb, label="Histogram")

  ## mark page as protected (i.e., cannot be closed)
  tag(win,"protected") <- append(tag(win,"protected"), length(nb))

  ## add functions for toolbar handlers.
  ## functions must be of type 'function(win)'
  pageno <- as.character(length(nb))
  tag(win,"run")[[pageno]] <- experim.ftable
  tag(win,"hist")[[pageno]] <- experim.hist.draw
  tag(win,"plot")[[pageno]] <- experim.hist.pvals.plot

  ## Distribution
  params.distr.gwl <- experim.hist.distr(main, win, ptests)

  ## we group frames for histogram and tests horizontally
  hist.grp <- ggroup(horizontal=TRUE, expand=TRUE, spacing=15, container=ptests)
  
  ## Histogram
  params.hist.gwl <- experim.hist.data(main, win, hist.grp)

  ## Goodness-of-Fit Tests
  params.gof.gwl <- experim.hist.gof(main, win, hist.grp)

  ## store list of widgets in window
  tag(win,"params.hist.gwl") <- c(params.distr.gwl, params.hist.gwl, params.gof.gwl)

  ## initialize total sample size
  experim.hist.update(win)

  ## return page
  return (ptests)
}


## --------------------------------------------------------------------------
## Histogram: Distribution 

experim.hist.distr <- function (main, win, group) {

  ## get data about generator object
  howdef <- tag(main,"distribution")$howdef

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## create a frame for inserting either the CDF or the quantile function
  frame <- gframe(" Distribution ", spacing=10, horizontal=FALSE, container=group)

  ## create a table for the layout of widgets
  tbl <- glayout(container=frame)
  tbl[1,2:3,anchor=c(-1,0)] <- "Need one of ..."

  ## construct character string for cdf and quantile function
  cdf.str <- icdf.str <- ""
  if (howdef == "built-in") { 
    ## get list of parameters (except domain) for built-in distribution
    distr.pm.gwl <- tag(main,"distr.pm.gwl")
    args <- paste(names(distr.pm.gwl),lapply(distr.pm.gwl,svalue), sep="=",collapse=", ")
    ## symbol for distribution object
    symbol <- tag(main,"distribution")$symbol
    ## R expression for evaluating cdf and inverse cdf
    if (exists(paste("p",symbol,sep=""), mode="function"))
      cdf.str  <- paste("p",symbol,"(x, ",args,")",sep="")                    
    if (exists(paste("q",symbol,sep=""), mode="function"))
      icdf.str <- paste("q",symbol,"(x, ",args,")",sep="")                    
  }

  ## create and place widgets
  cdf.cet  <- gcheckedit(label="CDF", checked=FALSE, text=cdf.str,
                         ttip="Cumulative distribution function of distribution",
                         coerce.with=wrap.function.body.NULL, width=50, container=tbl)

  tbl[2,2,anchor=c(-1,0)] <- params.gwl[["cdf.cbx"]] <- cdf.cet$cbx
  tbl[2,3,anchor=c(-1,0)] <- " <-  function(x) {"
  tbl[2,4,anchor=c(-1,0)] <- params.gwl[["cdf.edt"]] <- cdf.cet$edt
  tbl[2,5,anchor=c(-1,0)] <- "}"

  icdf.cet <- gcheckedit(label="quantile", checked=TRUE, text=icdf.str,
                         ttip="Qunatile function (Exact inverse cumulative distribution function)",
                         coerce.with=wrap.function.body.NULL, width=50, container=tbl)

  tbl[3,2,anchor=c(-1,0)] <- params.gwl[["icdf.cbx"]] <- icdf.cet$cbx
  tbl[3,3,anchor=c(-1,0)] <- " <-  function(x) {"
  tbl[3,4,anchor=c(-1,0)] <- params.gwl[["icdf.edt"]] <- icdf.cet$edt
  tbl[3,5,anchor=c(-1,0)] <- "}"

  ## (truncated) domain of distribution
  domain.gwl <- tag(main,"distr.params.gwl")
  domain.txt <- paste("truncated domain = (", svalue(domain.gwl[["lb"]]),",", svalue(domain.gwl[["ub"]]),")")
  params.gwl[["trunc.val"]] <- c(svalue(domain.gwl[["lb"]]), svalue(domain.gwl[["ub"]]))
  params.gwl[["istruncated"]] <- gcheckbox(domain.txt, checked=FALSE, container=tbl)
  tbl[4,2:4,anchor=c(-1,0)] <- params.gwl[["istruncated"]]
  tooltip(params.gwl[["istruncated"]]) <-
    paste("Check this box whenever CDF and/or quantile function have to be rescaled for the given domain.",
          "It is recommended to check 'CDF' in this case, too.", sep="\n")

  ## return list of widgets
  return (params.gwl)
}


## --------------------------------------------------------------------------
## Histogram: Data

experim.hist.data <- function (main, win, group) {
  
  ## store widgets for parameters in a list
  params.gwl <- list()

  ## create a frame for inserting parameters for histogram
  frame <- gframe(" Histogram ", spacing=10, horizontal=FALSE, container=group)

  ## create a table for the layout of widgets
  tbl <- glayout(container=frame)

  ## sample size for freqency table
  n.sample.edt   <- gedit("100000", coerce.with=as.integer, width=15, container=tbl)
  rep.sample.edt <- gedit("10",     coerce.with=as.integer, width=15, container=tbl)
  
  addHandlerKeystroke(n.sample.edt,   handler=function(h,...) {experim.hist.update(win)})
  addHandlerKeystroke(rep.sample.edt, handler=function(h,...) {experim.hist.update(win)})
  
  tbl[1,2,anchor=c(-1,0)] <- "sample size"
  tbl[1,3,anchor=c(-1,0)] <- params.gwl[["n.sample.edt"]] <- n.sample.edt
  tooltip(n.sample.edt) <- "Size of each random sample: a number between 1,000 and 10,000,000"
  
  tbl[2,2,anchor=c(-1,0)] <- "# of samples"
  tbl[2,3,anchor=c(-1,0)] <- params.gwl[["rep.sample.edt"]] <- rep.sample.edt 
  tooltip(rep.sample.edt) <- "Number of random samples to be drawn: a number between 10 and 1,000,000"
  
  ## number of break points for histogram
  breaks.edt <- gedit("101", coerce.with=as.integer, width=15, container=tbl)
  tbl[4,2,anchor=c(-1,0)] <- "# of break points"
  tbl[4,3,anchor=c(-1,0)] <- params.gwl[["breaks.edt"]] <- breaks.edt
  tooltip(breaks.edt) <- "Number of breakpoints for histogram (including boundary of domain): a number between 3 and 1001"

  ## return list of widgets
  return (params.gwl)

}

## --------------------------------------------------------------------------
## Histogram: Result of Goodness-of-Fit Tests

experim.hist.gof <- function(main, win, group) {

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## create a frame for displaying results of GoF tests
  frame <- gframe(" Goodness-of-Fit Tests ", spacing=10, horizontal=FALSE, container=group)
  
  ## create a table for the layout of widgets
  tbl <- glayout(container=frame)

  ## print total sample size
  tbl[1,2:3,anchor=c(-1,0)] <- "Total sample size:"
  tbl[1,4,anchor=c(-1,0)] <- params.gwl[["tot.sample.lbl"]] <- glabel("",container=tbl)

  ## print run time
  tbl[2,2:3,anchor=c(-1,0)] <- "Total run time:"
  tbl[2,4,anchor=c(-1,0)] <- params.gwl[["runtime.lbl"]] <- glabel("- - : - - . - - -",container=tbl)
  
  ## print p-values
  tbl[4,2:3,anchor=c(-1,0)] <- "p-value for ..."
  
  tbl[5,2,anchor=c(-1,0)] <- "[1]"
  tbl[5,3,anchor=c(-1,0)] <- "Chi-square test:"
  tbl[5,4,anchor=c(-1,0)] <- params.gwl[["pval.chisq.lbl"]] <- glabel("...", container=tbl)
  ## tbl[5,5,anchor=c(-1,0)] <- gimage("info",dirname="stock", container=tbl,
  ##           handler=show.help, action=list("rvgt.chisq") )

  tbl[6,2,anchor=c(-1,0)] <- "[2]"
  tbl[6,3,anchor=c(-1,0)] <- "M-Test:"
  tbl[6,4,anchor=c(-1,0)] <- params.gwl[["pval.Mtest.lbl"]] <- glabel("...", container=tbl)
  ## tbl[6,5,anchor=c(-1,0)] <- gimage("info",dirname="stock", container=tbl,
  ##           handler=show.help, action=list("rvgt.Mtest") )
  
  ## return list of widgets
  return (params.gwl)
}


## --------------------------------------------------------------------------
## Update histogram page

experim.hist.update <- function(win) {

  ## get data
  params.gwl <- tag(win,"params.hist.gwl")

  ## update total sample size
  val <- as.numeric(svalue(params.gwl[["n.sample.edt"]])) *
    as.numeric(svalue(params.gwl[["rep.sample.edt"]]))
  str <- if (isTRUE(val < 1e9)) prettyNum(as.integer(val), big.mark=" ") else prettyNum(val)
  svalue(params.gwl[["tot.sample.lbl"]]) <- str

  ## remove frequency table
  tag(win,"ftable") <- NULL
}
  

## --------------------------------------------------------------------------
## Compute RVGT frequency table and compute p-values for tests

experim.ftable <- function(win) {

  ## update statusbar
  ## sb <- tag(win,"sb")
  ## svalue(sb) <- "Computing histogram. This may take some time ..."

  ## list of widgets with data
  params.gwl <- tag(win,"params.hist.gwl")

  ## sample size
  n <- svalue(params.gwl[["n.sample.edt"]])
  rep <- svalue(params.gwl[["rep.sample.edt"]])
  
  if (n<1000 || n>1e7) {
    error.message("sample size < 1000 or > 1e7")
    tag(win,"ftable") <- NULL
    return(NULL)
  }
  if (rep<10 || rep>1e6) {
    error.message("# of samples < 10 or > 1e6")
    tag(win,"ftable") <- NULL
    return(NULL)
  }
  
  ## break points
  breaks <- svalue(params.gwl[["breaks.edt"]])

  if (breaks<3 || breaks>1001) {
    error.message("# of break points < 3 or > 1001")
    tag(win,"ftable") <- NULL
    return(NULL)
  }
 
  ## CDF and inverse CDF
  cdf.str  <- if( isTRUE(svalue(params.gwl[["cdf.cbx"]])))  svalue(params.gwl[["cdf.edt"]])  else NULL
  icdf.str <- if( isTRUE(svalue(params.gwl[["icdf.cbx"]]))) svalue(params.gwl[["icdf.edt"]]) else NULL

  if (is.null(cdf.str) && is.null(icdf.str)) {
    error.message("CDF / quantile function required")
    tag(win,"ftable") <- NULL
    return(NULL)
  }

  pdist <- tryCatch(eval(parse(text=paste(cdf.str))),
                    error=function(e){
                      error.message(e$message, title="UNU.RAN - Error: CDF"); NULL })

  qdist <- tryCatch(eval(parse(text=paste(icdf.str))),
                    error=function(e){
                      error.message(e$message, title="UNU.RAN - Error: quantile"); NULL })

  if (is.null(pdist) && is.null(qdist)) {
    error.message("Please select and/or insert 'CDF' or 'quantile'")
    tag(win,"ftable") <- NULL
    return (NULL);
  }

  ## truncated distribution
  if (isTRUE(svalue(params.gwl[["istruncated"]]))) {
    trunc <- params.gwl[["trunc.val"]]
  }
  else {
    trunc <- NULL
  }

  ## random variate generator
  rdist <- function(size) { ur(tag(win,"gen"), size) }

  ## start timer
  time.start <- proc.time()

  ## compute frequency table
  ftable <-
    tryCatch(rvgt.ftable(n=n, rep=rep, rdist=rdist, qdist=qdist, pdist=pdist,
                         trunc=trunc, breaks=breaks),
             error=function(e){
               error.message(e$message, title="UNU.RAN - Error: Cannot compute histogram")
               NULL
             } )

  if (is.null(ftable)) {
    tag(win,"ftable") <- NULL
    return(NULL)
  }
  
  ## stop timer
  run.time <- (proc.time() - time.start)[3]  ## "elapsed" time

  ## compute p-values for tests
  tag(win,"pvals.chisq") <- res.chisq <- rvgt.chisq(ftable)
  tag(win,"pvals.Mtest") <- res.Mtest <- rvgt.Mtest(ftable)
                        
  ## print result into frame
  svalue(params.gwl[["pval.chisq.lbl"]]) <- pretty.pval(res.chisq$pval)
  svalue(params.gwl[["pval.Mtest.lbl"]]) <- pretty.pval(res.Mtest$pval)

  ## print run time into frame
  min <- as.integer(floor(run.time / 60))
  sec <- as.integer(floor(run.time - min*60))
  msec <- as.integer(round(1000*(run.time - min*60 - sec)))
  svalue(params.gwl[["runtime.lbl"]]) <- sprintf("%02d:%02d.%03d", min, sec, msec)
  
  ## store ftable in window and reset counter for plots
  tag(win,"ftable") <- ftable
  tag(win,"plot.hist") <- FALSE
  tag(win,"plot.pval") <- FALSE
  
  ## update statusbar
  ## svalue(sb) <- "Histogram computed."

  ## return table
  return (ftable)
}                          
  
## --------------------------------------------------------------------------
## print p-value of last entry in list

pretty.pval <- function (pval) {
  p <- signif(pval[length(pval)],3)
  n.stars <- max( 0, abs(ceiling(log10(p))) - 1 )
  stars <- if (n.stars<=5) paste(rep("*",n.stars), collapse="") else "!!!!!"
  return (paste(p," ",stars))
}

## --------------------------------------------------------------------------
## draw histogram

experim.hist.draw <- function(win) {

  ## get ftable
  ftable <- tag(win,"ftable")
  if (is.null(ftable) || isTRUE(tag(win,"plot.hist"))) {
    ftable <- experim.ftable(win)
  }
  if (is.null(ftable)) return()
  
  ## draw histogram
  ggraphics(container=tag(win,"nb"), label="hist")
  Sys.sleep(0.5)
  plot(ftable)
  tag(win,"plot.hist") <- TRUE
}

## --------------------------------------------------------------------------
## plot p-values

experim.hist.pvals.plot <- function(win) {

  ## get ftable
  ftable <- tag(win,"ftable")
  if (is.null(ftable) || isTRUE(tag(win,"plot.pval"))) {
    ftable <- experim.ftable(win)
  }
  if (is.null(ftable)) return()

  ## list of p-values
  pvals <- list( tag(win,"pvals.chisq"), tag(win,"pvals.Mtest") )

  ## plot p-values
  ggraphics(container=tag(win,"nb"), label="p-values")
  Sys.sleep(0.5)
  plot.rvgt.htest(pvals, alpha = 0.001)
  tag(win,"plot.pval") <- TRUE
}

## --------------------------------------------------------------------------
