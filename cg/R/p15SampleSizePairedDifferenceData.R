## $Id:  $ 
## Paired Groups Case

## Sample Size Calculations for Paired Groups Data

setClass("cgPairedDifferenceSampleSizeTable",
         representation(ols.sstable="dataframeMatrixOrNULL",
                        settings="list"),
         prototype(ols.sstable=NULL,
                   settings=list()))

setMethod("samplesizeTable", "cgPairedDifferenceFit",
          samplesizeTable.cgPairedDifferenceFit <-
          function (fit, direction, mmdvec, power = 0.8, alpha = 0.05,
                    nmax = 1000, display = "print", ...) {
            ##
            ## PURPOSE:  Compute sample sizes based on variance
            ## estimation from Paired Difference Fit and other specifications
            ##
            ## NOTE: There is no capability to handle the Resistant & Robust
            ## fit since there is no decomposition of variance components
            ## analogous to the classical olsfit 
            ##
            ## Check input arguments
            dots <- list(...)
            validDotsArgs(dots, names="correction")
            correctionarg <- getDotsArgName(dots, "correction")
            if(!is.na(correctionarg)) {
              correction <- eval(parse(text=paste("dots$", correctionarg, sep="")))
              correction <- validArgMatch(correction, choices=c("df","none"))
            }
            else {
              correction <- "none"
            }
            
            direction <- validArgMatch(direction, c("increasing", "decreasing"))
            validNumeric(mmdvec, positive = TRUE)
            validAlpha(alpha)
            validPower(power)
            validNumeric(nmax, positive = TRUE, integer = TRUE)
            if (nmax < 2) {
              stop(cgMessage("The nmax argument needs to be an integer of 2 or greater."))
            }
            display <- validArgMatch(display, c("print","none","show"))

            ## The following argument defined in the generic is not needed,
            ## and thus ignored
            ngrpsarg <- getDotsArgName(dots, name="ngrps")
            if(!is.na(ngrpsarg)) {
              cat(cgMessage("The ngrps argument seems to be specified,",
                            "perhaps by an unnamed argument value. It is not applicable for this",
                            "method and is ignored...\n\n", warning=TRUE))
            }

            ## initializations
            settings <- fit@settings
            endptscale <- settings$endptscale
            planningname <- settings$analysisname
            endptname <- settings$endptname

            olsfit <- fit@olsfit

            ## Core calculations of inputs for noncentrality parameter (ncp)
            ## and denominator degrees of freedom (dendf)
            ## Note that the ngrps=2 argument is designed to be ignored in
            ## the body of the function.
            dendf <- function(n, ngrps=2) n - 1

            ## Calculate the noncentrality parameter
            if(correction=="none") {
              calc.ncp <- function(sigmaest, mmdvec, n) {
                ## variance components and relative efficiency index setup
                sigmaest.within <- sigmaest["within"]
                
                ## relative efficiency (in ratio scale, not percent scale)
                ## inverse of intraclass correlation
                releff <- (sigmaest["total"]/sigmaest.within)^2
                
                ncp <- (((n * (mmdvec)^2)/(2 * sigmaest.within^2)) *
                        releff)              
                ncp                                               
              }
            }
            else { ## if(correction=="df")
              calc.ncp <-  function(sigmaest, mmdvec, n) {

                fcndfcorrection <-  function(n) {
                  ## Fisher & Cochran / Cox degrees of freedom correction for using variance
                  ## estimates based on different degrees of freedom
                  ## See Fleiss (1986, pages 129 - 130)
                  ## The correction decreases the relative efficiency, and
                  ## and in turn the ncp as well. The sample size estimate
                  ## will increase.
                  ((2 * n + 1) * n)/((2 * n - 1) * (n + 2))
                }
                
                ## variance components and relative efficiency index setup
                sigmaest.within <- sigmaest["within"]
                
                ## relative efficiency (in ratio scale, not percent scale)
                ## inverse of intraclass correlation
                releff <- (sigmaest["total"]/sigmaest.within)^2

                ncp <- (((n * (mmdvec)^2)/(2 * sigmaest.within^2)) *
                        releff * fcndfcorrection(n))
                ncp                                               
              }
            }                       

            ols.sigmaest <- varianceTable(fit, display="none")@contents$"Spread(StdDev)"
            names(ols.sigmaest) <- c("within", "between", "total")

            ols.sstable <- samplesize(sigmaest=ols.sigmaest,
                                      endptscale=endptscale,
                                      planningname=planningname,
                                      endptname=endptname,
                                      ngrps=2, direction=direction,
                                      mmdvec=mmdvec,
                                      dendf=dendf,
                                      ncp=calc.ncp,
                                      alpha=alpha, power=power,
                                      nmax=nmax,
                                      display="none")

            settings <- c(settings, list(sigmaest = list(ols =
                                           ols.sigmaest),
                                         planningname = planningname,
                                         endptname=endptname,
                                         ngrps=2, direction=direction,
                                         mmdvec=mmdvec,
                                         dendf=dendf,
                                         alpha=alpha, power=power,
                                         nmax=nmax,
                                         display="none",
                                         correction=correction)
                          )

            returnObj <- new("cgPairedDifferenceSampleSizeTable", 
                             ols.sstable = ols.sstable,
                             settings = settings)
            
            if (display == "print") {
              print(returnObj)
            }
            else if (display == "show") {
              show(returnObj)
            }
            ## else display=="none"
            invisible(returnObj)
          }
          )

setMethod("print", "cgPairedDifferenceSampleSizeTable",
          print.cgPairedDifferenceSampleSizeTable <-
          function (x, title=NULL, endptname=NULL, ...)
          {
            ##
            ## PURPOSE: Semi-formatted print version of Sample Size Table
            ##
            ## NOTE: Had to use x as an argument because of the system defined
            ## generic. Would have preferred to use object; hence the first
            ## statement below.
            ##
            object <- x
            ## Input arguments check
            settings <- object@settings
            alpha <- settings$alpha
            power <- settings$power
            endptscale <- settings$endptscale
            ngrps <- 2
            nmax <- settings$nmax
            expunitname <- settings$expunitname
            correction <- settings$correction
            grpnames <- settings$grpnames
            
            sigmaest <- settings$sigmaest

            if (is.null(title)) {
              title <- paste("Sample Size Table from", settings$planningname)
            }
            else {
              validCharacter(title)
            }
            
            if (is.null(endptname)) {
              endptname <- settings$endptname
              if (!is.character(endptname)) {
                endptname <- ""
              }
            }
            else {
              validCharacter(endptname)
            }

            cat(title, "\n")
            if (endptname != "") {
              cat(paste("Endpoint:", endptname))
            }
            cat(paste(" (", grpnames[2], " vs. ", grpnames[1], ")", sep=""), "\n")
            
            ## Taken from base:::chartr help file
            .simpleCap <- function(x) {
              s <- strsplit(x, " ")[[1]]
              paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "",
                    collapse = " ")
            }

            diffmetric <- "Differences"
            if (settings$endptscale == "log") {
              diffmetric <- paste("Percent", .simpleCap(settings$direction),
                                  diffmetric)
            }
            else {
              diffmetric <- paste(.simpleCap(settings$direction), diffmetric)
            }
            cat(diffmetric, "\n")
            
            cat(paste(round(100 * power, 0), "% Power and ", round(100 *
                                                                   alpha, 0),
                      "% Significance Level\n", sep = ""))
            
            fmtdig <- function(x) {
              x$mmd <- fround(x$mmd, getNumDigits(x$mmd))
              rownames(x) <- x$mmd
              x$mmd <- NULL
              names(x) <- c(paste("n ", expunitname, "s", sep=""),
                            paste("N", "msmts", sep=" "))
              x[, 1] <- makeCensored(x[, 1], nmax)
              x[, 2] <- makeCensored(x[, 2], 2 * nmax)
              x
            }

            informSettings <- function(sigmaest) {
              cat("Variability Estimates ",  if (endptscale == "log") {
                "(Log scale)"},  "\n",
                  
                  paste("Within", expunitname), ": ",
                  signif(sigmaest["within"], 4), "\n",
                  
                  paste("Between", expunitname), ": ",
                  signif(sigmaest["between"], 4), "\n",
                  
                  paste("Total", expunitname), ": ",
                  signif(sigmaest["total"], 4), "\n\n",

                  sep="")
            }
            
            curwidth <- getOption("width")
            on.exit(options(width = curwidth), add=TRUE)
            if (curwidth < 500) {
              options(width = 500)
            }

            cat("\nClassical Least Squares Based",
                if(correction=="df") {
                  "(df Corrected)"
                }, "\n")
            informSettings(sigmaest$ols)
            print(fmtdig(as.data.frame(object@ols.sstable)), quote=FALSE)

            if(correction=="df") {
              cat("\nNOTE:")
              cat(paragraphWrap(paste("The df Corrected method is the\n",
                                      "Cochran / Cox degrees of freedom correction\n",
                                      "for using variances based on different degrees of freedom.",
                                      "\n\n")))
            }

            invisible()
          })

setMethod("show", "cgPairedDifferenceSampleSizeTable",
          show.cgPairedDifferenceSampleSizeTable <- function(object) showDefault(object))

setMethod("samplesizeGraph", "cgPairedDifferenceSampleSizeTable",
          samplesizeGraph.cgPairedDifferenceSampleSizeTable <-
          function(sstable, Nscale = "log", mmdscale = "log",
                   ...) {
            ##
            ## PURPOSE: create a graph of samplesize calculations
            ## for a Paired Difference design
            ##
            ## NOTE: Because of the special case of the estimate of the
            ## number of exp units, the existing
            ## samplesizegraph utility function is not used. So the following
            ## code is written to accomodate the axis labels and other text
            ## descriptions in the graph for the Paired Difference case.
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names = c("mmdticklabels", "nticklabels"))
            nscalearg <- getDotsArgName(dots, name="nscale")
            if(!is.na(nscalearg)) {
              nscale <- validArgMatch(nscale, c("log","original"))
              if(Nscale!=nscale) {
                cat(cgMessage("The nscale argument overrides the Nscale",
                              "argument if it is used and they differ", warning=TRUE))
              }
            }
            else {
              nscale <- Nscale
            }
            mmdscale <- validArgMatch(mmdscale, c("log","original"))

            ## The following arguments defined in the generic is not needed,
            ## and thus ignored
            devicearg <- getDotsArgName(dots, name="device")
            if(!is.na(devicearg)) {
              cat(cgMessage("The device argument seems to be specified,",
                            "perhaps by an unnamed argument value. It is not applicable for this",
                            "method and is ignored...\n\n", warning=TRUE))
            }
            cgthemearg <- getDotsArgName(dots, name="cgtheme")
            if(!is.na(cgthemearg)) {
              cat(cgMessage("The cgtheme argument seems to be specified,",
                            "perhaps by an unnamed argument value. It is not applicable for this",
                            "method and is ignored...\n\n", warning=TRUE))
            }

            settings <- sstable@settings
            endptscale <- settings$endptscale
            difftype <- if(endptscale == "log") "percent" else "simple"
            alpha <- settings$alpha
            power <- settings$power
            planningname <- settings$planningname
            endptname <- settings$endptname
            endptlabel <- makeEndptLabel(endptname, settings$endptunits)
            ngrps <- 2
            direction <- settings$direction
            nmax <- settings$nmax
            stamps <- settings$stamps
            sigmaest <- settings$sigmaest
            grpnames <- settings$grpnames
            correction <- settings$correction

            if(nmax < 2) {
              stop(cgMessage("The nmax argument needs to be an integer of 2 or greater."))
            }
            if(endptscale=="log") {
              if(direction=="decreasing" && any(mmdvec >= 100)) {
                stop(cgMessage("Percent Change Values of 100% or greater reduction",
                               "seem to have been entered.  Since this is",
                               "numerically impossible, please check and delete any such",
                               "values."))
              }
            }
            
            ols.sstable <- sstable@ols.sstable
            
            nticklabelsarg <- getDotsArgName(dots, "nticklabels")
            if (!is.na(nticklabelsarg)) {
              nticklabels <- eval(parse(text = paste("dots$", nticklabelsarg, sep = "")))
              validList(nticklabels, names = c("mod", "marks"),
                        argname = "nticklabels")
            }
            else {
              nticklabels <- NULL
            }
            mmdticklabelsarg <- getDotsArgName(dots, "mmdticklabels")
            if (!is.na(mmdticklabelsarg)) {
              mmdticklabels <- eval(parse(text = paste("dots$", mmdticklabelsarg, sep = "")))
              validList(mmdticklabels, names = c("mod", "marks"), argname = "mmdticklabels")
            }
            else {
              mmdticklabels <- NULL
            }

            if(length(mmdvec <- ols.sstable[, 1])==1) {
              stop(cgMessage("No graph is produced since only one",
                             "difference value was specified in the sample size",
                             "table"))
            }

            options(warn=-1)
            curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010)
            options(warn=0)
            on.exit(par(curpar))

            ## n refers to the number of experimental units for this method.
            ## N will be n * 2
            n <- ols.sstable[, 2]
            N <- ols.sstable[, 3]
            numberofgrps <- 2
            
            n.scaled <- scaleVar(n, endptscale=nscale)
            n.logscale <- if(nscale=="log") TRUE else FALSE
            
            mmd.logscale <- if(mmdscale=="log") TRUE else FALSE
            logandpct <- if(mmd.logscale &&  difftype=="percent") TRUE else FALSE
            mmd.scaled <- scaleVar(mmdvec, endptscale=mmdscale,
                                   percent=logandpct)

            parmar <- par(mar=c(5, 4, 4, 3) + 0.1)
            curpar$mar <- parmar$mar

            ## Graph the region line, but no axes content except the labels
            ## A key distinction is that the left side axis refers to the
            ## the number of experimental units, and the right side axis
            ## referes to the total number of measurements. 
            plot(mmd.scaled, n.scaled, type="n", pch=1, 
                 xlab="",
                 ylab=paste("Number of ",
                   settings$expunitname, "s", sep=""),
                 axes=FALSE,
                 xlim=rangeExtend(range(mmd.scaled),
                   pct=list(minside=6, maxside=0)))
            grid(lty=1)
            lines(mmd.scaled, n.scaled, type="b", pch=1)

            ## Now we work on the axes
            xlabchar <- paste(if(difftype=="percent") {
              ## "Percent Difference"
              if(direction=="decreasing") {
                "Percent REDUCTION"
              }
              else {
                "Percent INCREASE"
              }
            }
            else {
              ## "Simple Difference"
              if(direction=="decreasing") {
                "Simple REDUCTION"
              }
              else {
                "Simple INCREASE"
              }
            },
                              ": Detectable",
                              " Difference", sep="")

            mtext(side=1, text=xlabchar, line=2)
            if(is.expression(endptname)) {
              mtext(side=1, text=c(catCharExpr("in ", endptname),
                              paste(" (", grpnames[2], " vs. ", grpnames[1], ")", sep="")),
                              line=3)
            }
            else {
              mtext(side=1, text=paste("in", endptname,
                              paste(" (", grpnames[2], " vs. ", grpnames[1],
                                    ")", sep="") ),
                    line=3)
            }
            
            title(line = 2,
                  main = paste("Sample Size Graph\n", planningname, sep = ""),
                  cex = 1.1)
            mtext(side = 3, line = 1,
                  text = paste(round(100 * power, 0), " % Power ; ",
                    round(100 * alpha, 0), " % Significance Level ; ",
                    "Classical Variability Estimates ",
                    if (difftype == "percent") {
                      "(Log scale): "
                      },
                    sep = ""),
                  cex = 0.6)
            sigmaest.ols <- sigmaest$ols
            mtext(side = 3, line = 0.3,
                  text = paste("within=",
                    round(sigmaest$ols[1], 4),
                      "; between=",
                    round(sigmaest$ols[2], 4),
                    "; total=", round(sigmaest$ols[3], 4),
                    if(correction=="df") {
                      " (df correction used)"
                    },
                    sep = ""),
                  cex = 0.6)

            ## Filling in axes tick marks and such.
            ## First the y-axes
            names(n.scaled) <- as.character(n)
            n.tickmarks <- setupAxisTicks(n, logscale=n.logscale,
                                          digits=0)
            
            if(!is.null(nticklabels)) {
              validList(nticklabels, names=c("mod","marks"),
                        argname="nticklabels")
              n.tickmarks <- makeTickMarks(nticklabels, n.tickmarks)
            }
            
            N.ticklabels <- paste("", round(n.tickmarks*numberofgrps, 0))
            ycex <- yTicksCex(n.tickmarks, cexinit=0.8)

            axis(2, at=scaleVar(n.tickmarks, endptscale=nscale),
                   labels=names(n.tickmarks), adj=1, cex.axis=ycex, las=1)
            axis(4, at=scaleVar(n.tickmarks, endptscale=nscale),
                 labels=N.ticklabels, adj=0, cex.axis=min(ycex, 0.7), las=1)
            minmaxTicks(n, logscale=n.logscale, digits=0)
            mtext(side=4, text="Total Number of Measurements", line=1.8, cex=0.9, adj=0)
            
            mmdvec.ticks <- if(logandpct) {
              pctToRatio(mmdvec)
            }
            else {
              mmdvec
            }

            ## Now the x-axes
            mmdDigits <- getNumDigits(mmdvec)
            names(mmd.scaled) <- as.character(mmdvec)
            
            mmd.tickmarks <- setupAxisTicks(
                                            mmdvec.ticks,
                                            axis="x",
                                            logscale=mmd.logscale,
                                            digits=mmdDigits,
                                            ratio=logandpct,
                                            percent=logandpct,
                                            remticks=TRUE)
            
            if(length(mmd.tickmarks) < 3) {
              mmd.tickmarks <- setupAxisTicks(
                                              mmdvec.ticks,
                                              axis="x",
                                              logscale=mmd.logscale,
                                              digits=getNumDigits(mmdvec.ticks),
                                              ratio=logandpct,
                                              percent=logandpct,
                                              remticks=TRUE)
            }
  
            mmd.tickmarks <- mmd.tickmarks[names(mmd.tickmarks)!="0"]
  
            if(!is.null(mmdticklabels)) {
              validList(mmdticklabels, names=c("mod","marks"),
                        argname="mmdticklabels")
              mmd.tickmarks <- makeTickMarks(mmdticklabels, mmd.tickmarks,
                                             percent=logandpct)
            }

            xcex <- xTicksCex(mmd.tickmarks, cexinit=0.8)
            axis(1, at=scaleVar(mmd.tickmarks, endptscale=mmdscale),
                 labels=names(mmd.tickmarks), adj=1, cex.axis=xcex,
                 las=1)

            box()
            
              if(any(n >= nmax)) {
                .R. <- TRUE
                coordinates <- Hmisc::largest.empty(mmd.scaled, n.scaled,
                                                    width=0.30*(diff(range(
                                                      c(min(mmd.scaled),max(mmd.scaled))))),
                                                    height=0.25*(diff(range(
                                                      c(min(n.scaled),max(n.scaled))))),
                                                    method="area")
                legend(x=coordinates$rect$x[4], y=coordinates$rect$y[4], adj=0,
                       xjust=0, yjust=1, bty="n",
                       legend=paste("For at least one difference\n",
                         "the sample size\n",
                         "calculations were truncated\n",
                         "at",  paste(nmax, "per group."),
                         sep=" "),
                       text.col="red",
                       cex=0.70)
              }
                        
            usedgrid <- FALSE
            if (stamps) graphStampCG(grid = usedgrid)
            
            invisible()
          }
          )


