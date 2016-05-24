# $Id:  $
## Paired Groups Case

## Descriptive methods

## Experimental unit Profile graph
setMethod("profileGraph", "cgPairedDifferenceData",
          profileGraph.cgPairedDifferenceData <-
          function(data, ...) {
            ##
            ## PURPOSE: Graph paired experimental unit data points (profiles),
            ## using log scale
            ## annotations if log-scale analysis option is requested.
            ##
            ## Input argument handling
            dots <- list(...)
            validDotsArgs(dots, names=c("logscale", "ticklabels"))
            
            settings <- data@settings
            
            dfr.gcfmt <- data@dfr.gcfmt
            dfru <- data@dfru
            digits <- settings$digits
            ## potential override of logscale value
            ## from endptscale slot in cgPairedDifferenceData object
            logscalearg <- getDotsArgName(dots, "logscale")
            if(!is.na(logscalearg)) {
              logscale <- eval(parse(text=paste("dots$", logscalearg, sep="")))
              validBoolean(logscale)
            }
            else {
              logscale <- ifelse(settings$endptscale=="log", TRUE, FALSE)
            }

            ## End input argument handling
            offset <- settings$addconstant
            zeroscore <- settings$zeroscore
            endptlabel <- makeEndptLabel(settings$endptname, settings$endptunits)
            
            grpnames <- settings$grpnames
            
            options(warn=-1)
            if(logscale) {
              curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010,
                            mar=c(5, 4, 4, 3) + 0.1)
            }
            else {
              curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010)
            }
            options(warn=0)
            on.exit(par(curpar), add=TRUE)

            ## Convert the data to yet another form: a list of length = number of
            ## experimnetal units. Each component of the list will have length 2.
            endptlist <- as.data.frame(t(dfr.gcfmt[, 2:3]))
            names(endptlist) <- as.character(dfr.gcfmt[, 1])
            n <- length(endptlist)
            endpt <- unlist(endptlist)

            ## Processing needed for labelling lines later
            expunitlabels <- names(endptlist)
            
            grp1data <- sapply(endptlist, function(x) x[1])
            grp2data <- sapply(endptlist, function(x) x[2])
            grpdata <- list(grp1data, grp2data)

            ## Sorting is needed for staggering adjacent value labels
            grp1order <- order(grp1data)
            grp2order <- order(grp2data)

            ordexpunitlabelgrp1 <- expunitlabels[grp1order]
            ordexpunitlabelgrp2 <- expunitlabels[grp2order]

            stagexpunitlabelgrp1 <- stagexpunitlabelgrp2 <- vector("character",length=n)

            ## The key recipe:  For all even-indexed values, make sure the label extends
            ## out away from the ones before and after it.
            for(i in 1:n) {
              if(i%%2==0) {
                stagexpunitlabelgrp1[i] <-
                  paste(ordexpunitlabelgrp1[i],
                        paste(rep("-",
                                  max(nchar(c(ordexpunitlabelgrp1[i-1],
                                              ordexpunitlabelgrp1[i+1])),
                                      na.rm=T) * 2),
                              collapse=""),
                        "  ",
                        sep="")
                stagexpunitlabelgrp2[i] <-
                  paste("  ",
                        paste(rep("-",
                                  max(nchar(c(ordexpunitlabelgrp1[i-1],
                                              ordexpunitlabelgrp1[i+1])),
                                      na.rm=T) * 2),
                              collapse=""),
                        ordexpunitlabelgrp2[i],sep="")
              }
              else {
                stagexpunitlabelgrp1[i] <- paste(ordexpunitlabelgrp1[i], "  ", sep="")
                stagexpunitlabelgrp2[i] <- paste("  ", ordexpunitlabelgrp2[i], sep="")
              }
            }

            stagexpunitlabelgrp <- list(stagexpunitlabelgrp1, stagexpunitlabelgrp2)

            doProfile <- function(endptlist, endptlabel,
                                  grpdata, stagexpunitlabelgrp,
                                  analysisname, expunitname) {

              ## Set-up but no points plotted yet
              endpt <- unlist(endptlist)
              plot(1, median(endpt), xlim=c(0.5, 2.5),
                   ylim=c(min(endpt), max(endpt)), type="n", axes=FALSE,
                   xlab="", ylab=endptlabel,
                   main=paste(expunitname, " Profiles\n", analysisname, sep=""),
                   cex.main=1.1)

              ## Plot the actual data
              n <- length(endptlist)
              jitteredx <- jitter(rep(1:2, each=n))
              
              for(i in 1:n) {
                x <- jitteredx[c(i,(n+i))]
                y <- endptlist[[i]]
                lines(x=x, y=y, type="b", col=1)
              }

              ## Adding the expunit labels
              text(x=rep(min(jitteredx), n),
                   y=sort(grpdata[[1]]),
                   adj=1, labels=stagexpunitlabelgrp[[1]], col="red",
                   cex=0.6, xpd=TRUE)
              text(x=rep(max(jitteredx), n),
                   y=sort(grpdata[[2]]),
                   adj=0, labels=stagexpunitlabelgrp[[2]], col="red",
                   cex=0.6, xpd=TRUE)
              
              invisible()
            }
            
            doProfile(if(logscale) lapply(endptlist, log10) else endptlist,
                      endptlabel,
                      if(logscale) lapply(grpdata, log10) else grpdata,
                      stagexpunitlabelgrp,
                      settings$analysisname,
                      settings$expunitname)

            ## Create the y-axes
            tickmarks <- setupAxisTicks(endpt, logscale=logscale,
                                        digits=digits,
                                        offset=offset, ycex=0.8)
            ticklabelsarg <- getDotsArgName(dots, "ticklabels")
            if(!is.na(ticklabelsarg)) {
              ticklabels <- eval(parse(text=paste("dots$", ticklabelsarg,
                                         sep="")))
              tickmarks <- makeTickMarks(ticklabels, tickmarks, 
                                         offset=offset) 
            }
            
            if(logscale) { 
              mtext("log-spaced", side=2, line=2.25, cex=0.7)
              axis(2, at=log10(tickmarks), labels=names(tickmarks),
                   cex.axis=0.8, adj=1, las=1)

              log10tickmarks <- setupLog10AxisTicks(endpt, offset)              
              axis(4, at=log10tickmarks,
                   labels=names(log10tickmarks), cex.axis=0.7,
                   adj=0, las=1, tck=-0.0075)
              mtext(side=4, text=catCharExpr("Log10 scale of ", endptlabel),
                    line=2, cex=0.7, adj=0) 
            }

            else { ## if endpoint scale is original
              axis(2, at=tickmarks, labels=names(tickmarks),
                   cex.axis=0.8, adj=1, las=1)
            }

            ## Axes Customization
            grpnameticksettings <- setupGrpNameTicks(grpnames, 1:2)
            plotGrpNameTicks((grpnames), settings=grpnameticksettings)

            minmaxTicks(endpt, theaxis="y", logscale=logscale, digits=digits,
                        offset=offset, zeroscore=zeroscore)

            if(settings$stamps) graphStampCG(grid=FALSE)
            box()
            
            invisible()
          })



## Boxplot / Point Graph combination display of differences
setMethod("diffGraph", "cgPairedDifferenceData",
          diffGraph.cgPairedDifferenceData <- function(data, ...) {
            ##
            ## PURPOSE: Create a point graph / boxplot combination graph of the paired
            ## data experimental unit differences, using
            ## log scale annotations if log-scale analysis option is requested.
            ##
            settings <- data@settings
            
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names=c("logscale", "ticklabels"))

            ## potential override of logscale value
            ## from endptscale slot in cgOneFactorData object
            logscalearg <- getDotsArgName(dots, "logscale")
            if(!is.na(logscalearg)) {
              logscale <- eval(parse(text=paste("dots$", logscalearg, sep="")))
              validBoolean(logscale)
            }
            else {
              logscale <- ifelse(settings$endptscale=="log", TRUE, FALSE)
            }
            ##

            dfr.gcfmt <- data@dfr.gcfmt
            dfru <- data@dfru
            endptlabel <- makeEndptLabel(settings$endptname, settings$endptunits)

            grpf <- dfru$grpf
            endpt <- dfru$endpt
            grpnames <- settings$grpnames

            options(warn=-1)
            if(logscale) {
              curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010,
                            mar=c(5, 4, 4, 3) + 0.1)
            }
            else {
              curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010)
            }
            options(warn=0)
            on.exit(par(curpar), add=TRUE)

            diffendpt <- { if(logscale) dfr.gcfmt[, "difflogendpt"]
            else dfr.gcfmt[, "diffendpt"] }

            expunitlabels <- names(diffendpt) <- as.character(dfr.gcfmt[, 1])
            n <- length(expunitlabels)
            
            diffdataorder <- order(diffendpt)
            ordexpunitlabeldiffdata <- expunitlabels[diffdataorder]
            stagexpunitlabeldiffdata <- vector("character", length=n)

            ## The key recipe:  For all even-indexed values,
            ## make sure the label extends out
            ## away from the ones before and after it.
            for(i in 1:n) {
              if(i%%2==0) {
                stagexpunitlabeldiffdata[i] <-
                  paste(ordexpunitlabeldiffdata[i],
                        paste(rep(" ",
                                  max(nchar(c(ordexpunitlabeldiffdata[i-1],
                                              ordexpunitlabeldiffdata[i+1])),
                                      na.rm=T) * 2),
                              collapse=""),
                        "  ",
                        sep="")
              }
              else {
                stagexpunitlabeldiffdata[i] <-
                  paste(ordexpunitlabeldiffdata[i],
                        "  ", sep="")
              }
            }

            ## Helper Function
            doDiffPlot <- function(diffendpt, pars,
                                   stagexpunitlabel) {

              ## Set up the plot with boxplot and some dummy placeholders
              centers <- boxplot(NA, NA,
                                 diffendpt,
                                 pars=pars,
                                 axes=F, style.bxp="att",
                                 ylim = c(min(0, min(diffendpt)),
                                   max(diffendpt)))
              
              ## Reference Line at 0 = Equivalence
              abline(h=0, col=4, lty=1, lwd=2)
              
              ## Add mean to boxplot
              points(3, mean(diffendpt), pch=3)
              
              ## Add data
              points(jitter(rep(2, n)), diffendpt)

              ## Add experimental unit labels
              text(x=rep(1.25, n),
                   y=sort(diffendpt),
                   adj=1, labels=stagexpunitlabel, col=2, cex=0.8)

              invisible()
            }

            doDiffPlot(if(logscale) diffendpt/log(10) else diffendpt,
                       pars=list(medlty="blank", medpch=1, boxfill="NA",
                         yaxt="n"),
                       stagexpunitlabel=stagexpunitlabeldiffdata)

            refgrp <- settings$refgrp
            minuendgrp <- grpnames[grpnames != refgrp]

            ## Title
            title(main=paste("Individual Differences Distribution\n",
                    settings$analysisname,": ", minuendgrp, " vs. ",
                    refgrp, sep=""))

            ## Create the y-axes
            tickmarks <- setupAxisTicks(if(logscale) exp(diffendpt) else diffendpt,
                                        logscale=logscale,
                                        ycex=0.8,
                                        ratio=if(logscale) TRUE else FALSE,
                                        percent=if(logscale) TRUE else FALSE)
            ticklabelsarg <- getDotsArgName(dots, "ticklabels")
            if(!is.na(ticklabelsarg)) {
              ticklabels <- eval(parse(text=paste("dots$", ticklabelsarg,
                                         sep="")))
              tickmarks <- makeTickMarks(ticklabels, tickmarks) 
            }

            if(logscale) {
              mtext(paste("Percent Difference of", minuendgrp, "vs.",
                          refgrp), side=2, line=3)
              mtext("log-spaced", side=2, line=2.25, cex=0.7) 
              axis(2, at=log10(tickmarks),
                   labels=names(tickmarks),
                   cex.axis=0.8, adj=1, las=1)
              
              axis(4, at=log10(tickmarks),
                   labels=fmtRatio(tickmarks),
                   cex.axis=0.7,
                   adj=0, las=1, tck=-0.0075)
              mtext(side=4, text=paste("Ratio Difference of ",
                              minuendgrp, "vs.", refgrp),
                    line=2, cex=0.7, adj=0) 
            }
            else { ## if endpoint scale is original
              mtext(paste("Difference of", minuendgrp, "vs.",
                          refgrp), side=2, line=3)
              axis(2, at=tickmarks, labels=names(tickmarks),
                   cex.axis=0.8, adj=1, las=1)
            }

            ## Create the x-axis
            axis(1, 1:3, c(settings$expunitname, "Points", "Boxplot"))

            ## Last annotations
            minmaxTicks(if(logscale) exp(diffendpt) else diffendpt,
                        digits=settings$digits,
                        difference=FALSE,
                        theaxis="y", logscale=logscale,
                        percent=if(logscale) TRUE else FALSE)

            box()
            if(settings$stamps) graphStampCG(grid=FALSE)
            boxplotStamp()
            
            ## Add a warning if sample size is small
            if(n < 6) {
              mtext(side=1, line=3,
                    text=paste("*** NOTE: A boxplot is ineffective or useless when",
                      "the number of pairs is 5 or less ***"), col="red")
            }
            
            invisible()

          })


## Descriptive Table
setClass("cgPairedDifferenceDescriptiveTable",
         representation(contents="data.frame", settings="list"),
         prototype(contents=data.frame(), settings=list()))

setMethod("descriptiveTable", "cgPairedDifferenceData",
          descriptiveTable.cgPairedDifferenceData <- function(data, display="print", ...) {
            ##
            ## 
            ## PURPOSE: Create a table of quantiles, summary statistics of the data's
            ## 2 individual groups.
            ##
            ## NOTE: Currently this method is very similar to that for
            ## descriptiveTable.cgPairedDifferenceData. No censored data
            ## handling and the inclusion of a summary of differences are
            ## the main distinctions.
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names=c("logscale"))
            
            display <- validArgMatch(display, c("print","none","show"))

            settings <- data@settings
            ## potential override of logscale value
            ## from endptscale slot in cgOneFactorData object
            logscalearg <- getDotsArgName(dots, "logscale")
            if(!is.na(logscalearg)) {
              logscale <- eval(parse(text=paste("dots$", logscalearg, sep="")))
              validBoolean(logscale)
            }
            else {
              logscale <- ifelse(settings$endptscale=="log", TRUE, FALSE)
            }
            ## End input arguments check

            offset <- settings$addconstant
            zeroscore <- settings$zeroscore

            dlist <-  with(data@dfru, split(endpt, grpf))

            ## Add in differences
            dfr.gcfmt <- data@dfr.gcfmt
            if(logscale) {
              difflogendpt <- dfr.gcfmt[, "difflogendpt"] 
              diffendpt <- exp(difflogendpt)
              diffendpt.pct <- 100 * (diffendpt - 1)
              diffendpt.orig <- dfr.gcfmt[, "diffendpt"]
            }
            else {
              diffendpt <- dfr.gcfmt[, "diffendpt"]
            }

            sumsts <- function(x, sdandse=TRUE) {
              c(c(length(x),min(x),
                  quantile(x,c(0.25,0.50,0.75)),
                  max(x),mean(x)),
                if(sdandse) {
                  c(sd(x),stndErr(x))
                }
                else {
                  c(NA, NA)
                }
                )
            }

            comps <- sapply(dlist,
                            function(x) {
                              sumsts(x)
                            })
            if(!logscale) {
              comps <- cbind(comps, sumsts(diffendpt))
            }
            
            dimnames(comps)[[1]] <- c("n","Min","25%ile","Median","75%ile","Max",
                                      "Mean","StdDev","StdErr")
            if(logscale) {
              comps <- rbind(comps, sapply(dlist,
                                           function(x) {
                                             c(geoMean(x),
                                               geoMean(x)*stndErr(log(x)),
                                               stndErr(log(x))) } )
                             )
              selogvec <- comps[nrow(comps), ]
              if(any(selogvec > 0.50)) {
                warning(cgMessage("There is at least one",
                                  "group standard error in the log scale",
                                  "that exceeds 0.50, so the",
                                  "estimated geometric mean standard",
                                  "errors may be nonsensical and should be",
                                  "cautiously regarded.")
                        ) 
              }
              ## Drop that last row since it was only needed for a potential
              ## warning.
              comps <- comps[-nrow(comps), ]
              dimnames(comps)[[1]][(nrow(comps)-1):nrow(comps)] <- c("GeoMean",
                                                                     "SEGeoMean")
              
              if(!is.null(offset)) {
                for(i in c(2:7, 10)) {
                  comps[i, ] <- comps[i, ] - offset
                }
                mincol <- comps[2, ]
                mincol[mincol==min(mincol)] <- min(data@dfr)
                comps[2, ] <- mincol
              }
              
              if(!is.null(zeroscore)) {
                mincol <- comps[2, ]
                mincol[mincol==min(mincol)] <- 0
                comps[2, ] <- mincol
              }
              
              ## Add in summary statistics of differences
              ## Endpoints in Original, Ratio and Percent
              comps <- cbind(comps,
                             c(sumsts(diffendpt.orig, sdandse=TRUE), NA, NA),
                             c(sumsts(diffendpt, sdandse=FALSE),
                               if(logscale) geoMean(diffendpt) else NA,
                               if(logscale) geoMean(diffendpt)*stndErr(difflogendpt)),
                             c(sumsts(diffendpt.pct, sdandse=FALSE),
                               if(logscale) 100*(geoMean(diffendpt) - 1) else NA,
                               if(logscale)
                               100*geoMean(diffendpt)*stndErr(difflogendpt)
                               else NA)) 
            }

            comps <- as.data.frame(t(comps), stringsAsFactors=FALSE)
            grpnames <- settings$grpnames
            ordgrpnames <- orderPairedGroups(grpnames, settings$refgrp)
            
            if(logscale) {
              difflabels <- c(paste("Simple Difference:", paste(ordgrpnames[2], "-",
                                                    ordgrpnames[1], sep="")),
                              paste("Ratio Difference:", paste(ordgrpnames[2], "/",
                                                    ordgrpnames[1], sep="")),
                              paste("Percent Difference:",
                                    paste(ordgrpnames[2], " vs. ",
                                                                 ordgrpnames[1],
                                                                 sep="")))
              row.names(comps) <- c(grpnames, difflabels)
            }
            else {
              label <- paste("Difference:", paste(ordgrpnames[2], "-",
                                                  ordgrpnames[1], sep=""))
              row.names(comps) <- c(grpnames, label)
            }
            
            returnObj <- new("cgPairedDifferenceDescriptiveTable",
                             contents=comps,
                             settings=settings)
            
            if(display=="print") {
              print(returnObj)
            }
            else if(display=="show") {
              show(returnObj)
            }
            ## else display=="none"           
            invisible(returnObj)  
          } )


setMethod("print", "cgPairedDifferenceDescriptiveTable",
          print.cgPairedDifferenceDescriptiveTable <-
          function(x, digits=NULL, title=NULL, endptname=NULL, ...) {
            ##
            ## PURPOSE: Semi-formatted print version of Descriptive Table
            ## 
            ## NOTE: Had to use x as an argument because of the system defined
            ## generic. Would have preferred to use object; hence the first
            ## statement below.
            ##
            object <- x
            settings <- object@settings

            if(is.null(digits)) {
              digits <- settings$digits
            }
            else {
              digits <- validArgDigits(digits)
            }
            curscipen <- getOption("scipen")
            on.exit(options(scipen=curscipen), add=TRUE)
            options(scipen=9)            
            
            if(is.null(title)) {
              title <- paste("Descriptive Table of", settings$analysisname) 
            }
            else {
              validCharacter(title)
            }
            if(is.null(endptname)) {
              endptname <- settings$endptname
              if(!is.character(endptname)) {
                endptname <- ""
              }
            }
            else {
              validCharacter(endptname)
            }

            cat(title,"\n")
            if(endptname!="") { cat(paste("Endpoint:", endptname, "\n")) }

            out <- object@contents
            fmt.out <- out
            
            if(nrow(fmt.out)==5) { ## log scale, ratio and percent
              fmt.out[1:3, -1] <- fround(out[1:3, -1], digits=digits)
              fmt.out[4, c(2:7, 10:11)] <- fmtRatio(unwind(out[4, c(2:7,
                                                                    10:11)]))
              fmt.out[5, c(2:7, 10:11)] <- fmtPercent(unwind(out[5, c(2:7,
                                                                      10:11)]))
              fmt.out[4:5, 8:9] <- "NA"
            }
            else { ## original scale
              fmt.out[, -1] <- fround(out[, -1], digits=digits)
            }
            
            fmt.out[,"n"] <- fround(as.numeric(out[, "n"]), 0)
            
            curwidth <- getOption("width")
            on.exit(options(width=curwidth), add=TRUE)
            if(curwidth < 500) { options(width=500)}
            print(fmt.out, quote=FALSE)
            invisible()
          })

setMethod("show", "cgPairedDifferenceDescriptiveTable",
          show.cgPairedDifferenceDescriptiveTable <- function(object) showDefault(object))


## Correlation Table

setClass("cgPairedDifferenceCorrelationTable",
         representation(contents="data.frame", settings="list"),
         prototype(contents=data.frame(), settings=list()))

setMethod("correlationTable", "cgPairedDifferenceData",
          correlationTable.cgPairedDifferenceData <- function(data,
                                                              display="print", ...) {
            ##
            ## PURPOSE: Create a table of the estimated correlation between the pairs
            ##
            ## Compute Pearson and Spearman indices, with log scale if option
            ## is chosen
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names=c("logscale"))
            
            display <- validArgMatch(display, c("print","none","show"))

            settings <- data@settings
            ## potential override of logscale value
            ## from endptscale slot in cgPairedDifferenceData object
            logscalearg <- getDotsArgName(dots, "logscale")
            if(!is.na(logscalearg)) {
              logscale <- eval(parse(text=paste("dots$", logscalearg, sep="")))
              validBoolean(logscale)
            }
            else {
              logscale <- ifelse(settings$endptscale=="log", TRUE, FALSE)
            }
            ## End input arguments check
            dfr.gcfmt <- data@dfr.gcfmt
            grp1 <- dfr.gcfmt[, 2]
            grp2 <- dfr.gcfmt[, 3]
            names(grp1) <- names(grp2) <- dfr.gcfmt[, 1]
            
            ## Original Scale Results which are always computed
            pearson.orig <- cor.test(grp1, grp2, method="pearson")
            
            ## For Spearman, scale does not matter, since log is a monotone
            ## function
            
            spearman.results <- suppressWarnings(cor.test(grp1, grp2, method="spearman"))
            
            ## Table construction
            if(logscale) {
              pearson.log <- cor.test(log(grp1), log(grp2), method="pearson")
              comps <- data.frame(correlation=
                                  c(pearson.orig$estimate,
                                    pearson.log$estimate,
                                    spearman.results$estimate))
              row.names(comps) <- c("Pearson Original", "Pearson Log",
                                    "Spearman")
            }
            
            else { ## original scale
              comps <- data.frame(correlation=
                                  c(pearson.orig$estimate,
                                    spearman.results$estimate))
              row.names(comps) <- c("Pearson","Spearman")
            }
            
            returnObj <- new("cgPairedDifferenceCorrelationTable",
                             contents=comps,
                             settings=settings)
            
            if(display=="print") {
              print(returnObj)  
            }
            else if(display=="show") {
              show(returnObj)
            }
            ## else display=="none"                
            invisible(returnObj)  
          })


setMethod("print", "cgPairedDifferenceCorrelationTable",
          print.cgPairedDifferenceCorrelationTable <-
          function(x, title=NULL, endptname=NULL, ...) {
            ##
            ## PURPOSE: Semi-formatted print version of Correlation Table
            ## 
            ## NOTE: Had to use x as an argument because of the system defined
            ## generic. Would have preferred to use object; hence the first
            ## statement below.
            ##
            object <- x
            settings <- object@settings
            
            if(is.null(title)) {
              title <- paste("Correlation Table of", settings$analysisname) 
            }
            else {
              validCharacter(title)
            }
            if(is.null(endptname)) {
              endptname <- settings$endptname
              if(!is.character(endptname)) {
                endptname <- ""
              }
            }
            else {
              validCharacter(endptname)
            }

            cat(title,"\n")
            cat(paste(settings$grpnames, collapse=" and "), "\n")
            if(endptname!="") { cat(paste("Endpoint:", endptname, "\n\n")) }
            
            fmt.out <- object@contents
            fmt.out[, "correlation"] <- fround(fmt.out[, "correlation"], 2)
            
            curwidth <- getOption("width")
            on.exit(options(width=curwidth), add=TRUE)
            if(curwidth < 500) { options(width=500)}
            print(fmt.out, quote=FALSE)
            invisible()
          })

setMethod("show", "cgPairedDifferenceCorrelationTable",
          show.cgPairedDifferenceCorrelationTable <- function(object) showDefault(object))
