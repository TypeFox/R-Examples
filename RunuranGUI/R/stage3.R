#############################################################################
##
##  Stage 3: Notebook for making some tests with UNU.RAN generator
##
#############################################################################
##
##  Synopsis:
##
##    stage3 ( main, gen )
##
##  Arguments:
##    main ... main window. It must contain all information created in
##             stages 1 and 2.
##    gen  ... UNU.RAN generator object
##
##
##  All data about distribution are stored in 'main'.
##  All data about this notebook are stored in 'win'.
##
##  Additional pages:
##    hist.R   ... histogram based tests
##    aerror.R ... approximation errors for numerical inversion
##
#############################################################################

## --------------------------------------------------------------------------

stage3 <- function(main, gen) {

  ## get data about generator object
  distr.str  <- tag(main,"distribution")$name
  method.str <- tag(main,"method")$symbol

  ## create new window with notebook
  win <- gwindow(paste("Runuran:",distr.str,"-",method.str))
  nb <- experim.notebook(win)

  ## store generator object in window
  tag(win,"gen") <- gen
  
  ## create page with R code
  experim.Rcode(main, win, nb)

  ## create page with properties of generator object
  experim.properties(main, win, nb)

  ## create page for histogram based tests
  experim.histtest(main, win, nb)

  ## create page for computing approximation errors of numerical inversion
  experim.inverror(main, win, nb)

  ## add statusbar
  ## tag(win,"sb") <- gstatusbar("", container=win)
  
  ## show first page
  svalue(nb) <- 1
  
}


## --------------------------------------------------------------------------
## Create notebook with toolbar

experim.notebook <- function (win) {

  group <- ggroup(horizontal=FALSE, spacing=10, container=win)

  ## toolbar buttons:

  ## quit notebook
  aQuit  <- gaction(label="Quit",   icon="quit",
                    tooltip="Quit notebook",
                    handler=function(h,...){dispose(win)})

  ## close current page
  aClose <- gaction(label="Close", icon="close",
                    tooltip="Close current page",
                    handler=function(h,...) {
                      nb <- tag(win,"nb")
                      pageno <- tag(win,"pageno")
                      if(! (pageno %in% tag(win,"protected"))) dispose(nb)
                    })
  
  ## evaluate page
  aRun <- gaction(label="Run", icon="execute",
                  tooltip="Evaluate page",
                  handler=function(h,...) {
                    nb <- tag(win,"nb")
                    pageno <- tag(win,"pageno")
                    func <- tag(win,"run")[[pageno]]
                    if(!is.null(func)) func(win)
                  })
  enabled(aRun) <- FALSE

  ## draw histogram
  aHist <- gaction(label="Histogram", icon="hist",
                  tooltip="Draw histogram",
                  handler=function(h,...) {
                    nb <- tag(win,"nb")
                    pageno <- tag(win,"pageno")
                    func <- tag(win,"hist")[[pageno]]
                    if(!is.null(func)) func(win)
                  })
  enabled(aHist) <- FALSE

  ## plot result
  aPlot <- gaction(label="Plot", icon="plot",
                  tooltip="Plot results",
                  handler=function(h,...) {
                    nb <- tag(win,"nb")
                    pageno <- tag(win,"pageno")
                    func <- tag(win,"plot")[[pageno]]
                    if(!is.null(func)) func(win)
                  })
  enabled(aPlot) <- FALSE

  
  ## create toolbar
  tb <- gtoolbar(list(Quit=aQuit,Close=aClose,Run=aRun,Hist=aHist,Plot=aPlot),
                 container=group)

  ## create notebook
  nb <- gnotebook(container=group, expand=TRUE, closebuttons = FALSE)

  ## store notebook in window
  tag(win,"nb") <- nb
  
  ## protected pages
  ## (we protect at least the first page)
  tag(win,"protected") <- 1

  ## list of handlers for buttons "Run", "Histogram" and "Plot"
  ## list entries must be functions of type 'function(win)'
  tag(win,"run")  <- list()
  tag(win,"hist") <- list()
  tag(win,"plot") <- list()

  ## store initial page number
  tag(win,"pageno") <- as.character(length(nb))

  ## update toolbar
  addHandlerChanged(nb, handler=function(h,...) {
    ## enable or disable buttons
    pageno <- ifelse(is.null(h$pageno), svalue(h$obj), h$pageno)
    pageno <- as.character(pageno)
    enabled(aClose) <- if(pageno %in% tag(win,"protected"))  FALSE else TRUE
    enabled(aRun)  <- if(is.null(tag(win,"run")[[pageno]]))  FALSE else TRUE
    enabled(aHist) <- if(is.null(tag(win,"hist")[[pageno]])) FALSE else TRUE
    enabled(aPlot) <- if(is.null(tag(win,"plot")[[pageno]])) FALSE else TRUE
    ## store current page number
    tag(win,"pageno") <- pageno
  } )

  ## return handler for notebook
  return(nb)
}
  

## --------------------------------------------------------------------------
## Create page with R code

experim.Rcode <- function (main, win, nb) {

  ## compose R code
  Rcode <- param.get.Rcode(main, internal=FALSE)

  ## compile text
  text <- paste("\n",
                "## Create UNU.RAN distribution object\n", Rcode$distr, "\n\n",
                "## Create UNU.RAN generator object\n", Rcode$gen, "\n\n",
                "## Draw sample\n", Rcode$sample, "\n\n",
                sep="")

  ## add page for R code
  pcode <- gtext(text=text, label="R Code", container=nb)
  tag(win,"protected") <- append(tag(win,"protected"), length(nb))
  
  ## return page
  return (pcode)
}


## --------------------------------------------------------------------------
## Create page with properties of generator object

experim.properties <- function (main, win, nb) {

  ## get string with properties
  ## prop <- capture.output(unuran.details( tag(win,"gen") ))
  prop <- capture.output(show( tag(win,"gen") ))

  ## compile text
  text <- paste(prop,sep="\n")

  ## add page for R code
  pprop <- gtext(text=text, label="Properties", container=nb)
  tag(win,"protected") <- append(tag(win,"protected"), length(nb))
  
  ## return page
  return (pprop)
}

## --------------------------------------------------------------------------
