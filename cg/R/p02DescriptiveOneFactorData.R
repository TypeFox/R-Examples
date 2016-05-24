## $Id: p02DescriptiveOneFactorData.R 4891 2014-02-21 17:49:45Z bpikouni $
## One-Factor Unpaired Groups Case

## Descriptive methods for One-Factor Unpaired Groups Data

## Point Graph
setMethod("pointGraph", "cgOneFactorData",
          pointGraph.cgOneFactorData <-
          function(data, ...) {
            ##
            ## PURPOSE: Create a point graph of the data, using log scale annotations
            ## if log-scale analysis option is requested. What I call "Point
            ## Graph" is also known as a dotplot, stripplot, or scattergram.
            ##
            ## Input argument handling
            dots <- list(...)
            validDotsArgs(dots, names=c("logscale", "ticklabels"))
            
            settings <- data@settings
            
            dfru <- data@dfru
            digits <- settings$digits
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
            ## End input argument handling

            offset <- settings$addconstant
            zeroscore <- settings$zeroscore
            endptlabel <- makeEndptLabel(settings$endptname, settings$endptunits)
            
            grpf <- dfru$grpf
            endpt <- dfru$endpt
            grpnames <- settings$grpnames
            numberofgrps <- length(grpnames)

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

            has.censored <- data@has.censored
            if(has.censored) {
              status <- dfru$status
            }

            doPlot <- function(x, y, endptlabel,
                               numberofgrps,
                               analysisname, has.censored) {
              plot((as.numeric(x)), y,
                   ylab=endptlabel,  xlab="",
                   xlim=c(0, numberofgrps) + 0.5, axes=FALSE,
                   main=paste("Point Graph\n", analysisname, sep=""),
                   cex.main=1.1, type="n")
              if(has.censored) {
                xn <- jitter(as.numeric(x))
                if(sum(status==0) > 0) {
                  text((xn)[status==0], y[status==0], labels=">", srt=90,
                       col="black")
                }
                points((xn)[status==1], y[status==1])
                if(sum(status==2) > 0) {
                  text((xn)[status==2], y[status==2], labels="<", srt=90,
                       col="black")
                }
                if(any(status==3)) {
                  stop(cgMessage("There is at least one",
                                 "interval censored value",
                                 "and no pointGraph method",
                                 "is currently available to handle these.")) 
                }
              }
              else {
                points(jitter(as.numeric(x)), y)
              }
              invisible()
            }

            doPlot(grpf, if(logscale) log10(endpt) else endpt,
                   endptlabel,
                   numberofgrps,
                   settings$analysisname,
                   has.censored)

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
            grpnameticksettings <- setupGrpNameTicks(grpnames, 1:numberofgrps)
            plotGrpNameTicks((grpnames), settings=grpnameticksettings)

            minmaxTicks(endpt, theaxis="y", logscale=logscale, digits=digits,
                        offset=offset, zeroscore=zeroscore)

            if(settings$stamps) graphStampCG(grid=FALSE)
            box()

            ## faint gray vertical dividers between groups
            abline(v=seq(1.5, numberofgrps-0.5, 1), col="lightgray")
            
            invisible()
          }
          )

## Boxplot
setMethod("boxplot", "cgOneFactorData",
          boxplot.cgOneFactorData <-
          function(x, ...) {
            ##
            ## PURPOSE: Create a boxplot graph of the data, using log scale annotations
            ## if log-scale analysis option is requested.
            ##
            ## NOTE: Had to use x as an argument because of the system defined
            ## generic. Would have preferred to use data; hence the first
            ## assignment below.
            data <- x
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

            dfru <- data@dfru
            digits <- settings$digits
            offset <- settings$addconstant
            zeroscore <- settings$zeroscore
            endptlabel <- makeEndptLabel(settings$endptname, settings$endptunits)

            grpf <- dfru$grpf
            endpt <- dfru$endpt
            grpnames <- settings$grpnames
            numberofgrps <- length(grpnames)

            options(warn=-1)
            if(logscale) {
              curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010,
                            mar=c(5, 4, 4, 3) + 0.1)
            }
            else {
              curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010)
            }
            options(warn=0)
            on.exit(par(curpar))

            has.censored <- data@has.censored 

            bppars <- list(medlty="blank", medpch=1, boxfill="NA",
                           yaxt="n")

            ## Helper Function
            doBoxPlot <- function(dfru, logscale, pars, endptlabel, numberofgrps,
                                  analysisname, has.censored) {
              if(!has.censored) {
                bxplist <- with(dfru, split(endpt, grpf))
                if(logscale) log10bxplist <- lapply(bxplist, log10)
                
                boxplot(if(logscale) log10bxplist else bxplist,
                        pars=pars,
                        ylab=endptlabel,
                        names=rep("", numberofgrps),
                        main=paste("Boxplot Graph\n",
                          analysisname, sep=""), cex.main=1.1)
                ## add means
                points(1:numberofgrps,
                       unlist(lapply(if(logscale) { log10bxplist } else {bxplist},
                                     function(x) {
                                       mean(x, na.rm=TRUE)
                                     }
                                     )), pch=3)
                ## Print some warnings if sample sizes are too small
                n <- sapply(bxplist, length)
                if(all(n < 6)) {
                  mtext(side=1, line=3,
                        text=paste("*** NOTE: Boxplots are ineffective when",
                          "all sample sizes are 5 or less ***"), col="red", cex=0.7)
                }
                
                else {
                  if(any(n < 6)) {
                    mtext(side=1, line=3,
                          text=paste("*** NOTE: At least one group has a sample size",
                            "of 5 or less,\nso its boxplot representation should be",
                            "ignored."), col="red", cex=0.7)
                  }
                }
              }
              else if(has.censored) {
                boxplotcensoreddata(dfru,
                                    logscale=logscale,
                                    pars=pars,
                                    ylab=endptlabel,
                                    names=rep("", numberofgrps),
                                    main=paste("Boxplot Graph\n",
                                      analysisname, sep=""), cex.main=1.1)
                ## no means added
              }
              invisible()
            }

            doBoxPlot(dfru, logscale, bppars, endptlabel, numberofgrps,
                      settings$analysisname, has.censored)

            ## Create the y-axes
            tickmarks <- setupAxisTicks(endpt, logscale=logscale,
                                        digits=digits,
                                        offset=offset,
                                        ycex=0.8)
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
            grpnameticksettings <- setupGrpNameTicks(grpnames, 1:numberofgrps)
            plotGrpNameTicks((grpnames), settings=grpnameticksettings)

            minmaxTicks(endpt, theaxis="y", logscale=logscale, digits=digits,
                        offset=offset, zeroscore=zeroscore)

            box()
            if(settings$stamps) graphStampCG(grid=FALSE)
            if(!has.censored) boxplotStamp()
            
            invisible()
          })


boxplotcensoreddata <- function(dfru, logscale, pars, ylab,
                                names, main, ...) {
  ## PURPOSE: Experimental Function for boxplot graph when censored data
  ## is present.
  ##
  ## NOTE: Currently this is designed only for use within the boxplot
  ## method for one factor data. Therefore minimal input checking is done.

  ## status <- dfru$status
  ## Draw off the ideas proposed by Gentleman & Crowley (1991)
  ## use KM type estimates for quantiles

  ## First need to construct the "boxplot.stats" structure used
  ## boxplot.default via bxp
  options(warn=-1)
  scfit <-  survfit(survival::Surv(endpt1, endpt2,
                                   status, "interval") ~ -1 + grpf,
                    data=dfru)
  summ.scfit <- summary(scfit, censored=TRUE)
  options(warn=0)

  ## Use estimated survival distribution function
  ## values to get quantiles
  scfit.dfr <- with(summ.scfit,
                    ## The factors in summ.scfit$grpf sometimes have "grpf="
                    ## prepended to the original factors from dfru$grpf
                    ## Remove that (if it is there) and 
                    ## refactor them using the original factor levels.
                    data.frame(grpf=factor(sub ("^grpf=", "", strata),
                                 levels=levels(dfru$grpf)),
                               endpt=time,
                               sdf.prob=surv,
                               cdf.prob=1-surv))
  
  scfit.dfr <- merge(dfru, scfit.dfr, by=c("grpf","endpt"))
  scfit.dfr <- scfit.dfr[order(scfit.dfr$grpf), ]
  grpnames <-  levels(dfru$grpf)

  bx.p <- lapply(split(scfit.dfr, scfit.dfr$grpf),
                 function(x) {
                   y <- x[x$status==1, ]
                   
                   thegrpf <- as.character(unique(x$grpf))
                   ## Follow survfit conventions for estimating
                   ## the median. 
                   ## 
                   p50 <- qminmin(x$sdf.prob, x$endpt, 0.50)
                   p25 <- qminmin(x$sdf.prob, x$endpt, 0.25)
                   p75 <- qminmin(x$sdf.prob, x$endpt, 0.75)

                   ## Initialize the whisker extremes
                   ## with the usual conventions
                   ## (see stats:::fivenum)

                   ## uncensored and left censored on lower end
                   lwe <- min(x$endpt[is.element(x$status, c(1,2))])

                   ## uncensored and right censored on upper end
                   uwe <- max(x$endpt[is.element(x$status, c(1,0))]) 

                   ## since there is censoring, then we may have
                   ## NA's
                   bxpstats <- c(lwe, p25, p50, p75, uwe)
                   nas <- is.na(bxpstats)
                   
                   if(nas[3]) {
                     warning(cgMessage("No boxplot is possible for the",
                                       thegrpf, "group",
                                       "since the median could",
                                       "not be calculated due to",
                                       "the amount of",
                                       "censoring or lack of",
                                       "distinct complete observations.",
                                       warning=TRUE))
                     return(list(stats=rep(NA, 5)))
                   }

                   if(nas[2]) {
                     warning(cgMessage("The 25th Percentile",
                                       "for the", 
                                       thegrpf, "group",
                                       "boxplot", 
                                       "could not be calculated due to",
                                       "the amount of",
                                       "censoring or lack of",
                                       "distinct complete observations.",
                                       warning=TRUE))
                     bxpstats[2] <- p50 ## set for plot purposes
                     ## and set the bottom whisker to the
                     ## the 'lowest' (and possibly left censored) value
                     lower.ext <- bxpstats[1] <- lwe
                   }
                   else {
                     lower.ext <- p25 - 1.5 * (p75 - p25)                   
                   }
                   
                   if(nas[4]) {
                     warning(cgMessage("The 75th Percentile",
                                       "for the", 
                                       thegrpf, "group",
                                       "boxplot", 
                                       "could not be calculated due to",
                                       "the amount of",
                                       "censoring or lack of",
                                       "distinct complete observations.",
                                       warning=TRUE))
                     bxpstats[4] <- p50 ## set for plot purposes
                     ## and set the top whisker to the
                     ## the 'highest' (and possibly right censored) value
                     upper.ext <- bxpstats[5] <- uwe
                   }
                   else {
                     upper.ext <- p75 + 1.5 * (p75 - p25)
                   }

                   ## follow and annotate boxplot.stats convention
                   allendpt <- x$endpt
                   allstatus <- x$status

                   out.indx <- (allendpt < lower.ext | allendpt > upper.ext)
                   
                   if(sum(out.indx) > 0) {
                     out.endpt <- allendpt[out.indx]
                     out.status <- rep(allstatus[out.indx]) 
                     out.group <- rep(thegrpf, length(out.endpt))
                     out.groupnumber <- match(out.group, grpnames)
                     ## for any tied endpt values jitter
                     ## groupnumber for later graph separation
                     out.jitter <- duplicated(out.endpt)
                     if(sum(out.jitter) > 0) {
                       out.groupnumber[out.jitter] <-
                         jitter(out.groupnumber[out.jitter])
                     }
                     bxpstats[1] <- min(allendpt[!(out.indx &
                                                   allstatus!=0)])
                     bxpstats[5] <- max(allendpt[!(out.indx &
                                                   allstatus!=2)])
                   }
                   else {
                     out.endpt <- out.status <- out.group <-
                       out.groupnumber <- NULL
                   }
                   
                   ## at box edges (25th and 75th percentiles)
                   ## or whiskers
                   allcensendpt <- allendpt[allstatus!=1]
                   allcensstatus <- allstatus[allstatus!=1]
                   
                   edges.indx <- (is.element(allcensendpt,
                                             c(p25, p75,
                                               bxpstats[c(1,5)])))
                   if(sum(edges.indx) > 0) {
                     edges.endpt <- allcensendpt[edges.indx]
                     edges.status <- rep(allcensstatus[edges.indx]) 
                     edges.group <- rep(thegrpf, length(edges.endpt))
                     edges.groupnumber <- match(edges.group,
                                                grpnames)
                     ## for any tied endpt values jitter
                     ## groupnumber for later graph separation
                     edges.jitter <- duplicated(edges.endpt)
                     if(sum(edges.jitter) > 0) {
                       edges.groupnumber[edges.jitter] <-
                         jitter(edges.groupnumber[edges.jitter])
                     }
                   }
                   else {
                     edges.endpt <- edges.status <- edges.group <-
                       edges.groupnumber <- NULL
                   }

                   return(list(stats=bxpstats,
                               n=length(allendpt),
                               out.endpt=out.endpt, out.status=out.status,
                               out.group=out.group,
                               out.groupnumber=out.groupnumber,
                               edges.endpt=edges.endpt,
                               edges.status=edges.status,
                               edges.group=edges.group,
                               edges.groupnumber=edges.groupnumber))
                 })

  ## adapted construction since any individual unusual points
  ## will be added separately after the bxp() graph call.
  bxp.sts <- list(stats=sapply(bx.p, function(x) x$stats),
                  n=sapply(bx.p, function(x) x$n),
                  conf=numeric(0L),
                  out=numeric(0L),
                  group=numeric(0L),
                  names=rep("", length(grpnames)),
                  out.endpt=unlist(sapply(bx.p, function(x) x$out.endpt)),
                  out.status=unlist(sapply(bx.p,
                    function(x) x$ out.status)),
                  out.group=unlist(sapply(bx.p, function(x)
                    x$out.group)),
                  out.groupnumber=unlist(sapply(bx.p, function(x)
                    x$out.groupnumber)),
                  edges.endpt=unlist(sapply(bx.p,
                    function(x) x$edges.endpt)),
                  edges.status=unlist(sapply(bx.p,
                    function(x) x$ edges.status)),
                  edges.group=unlist(sapply(bx.p, function(x)
                    x$edges.group)),
                  edges.groupnumber=unlist(sapply(bx.p, function(x)
                    x$edges.groupnumber))
                  )
  

  if(logscale) { bxp.sts$stats <- log10(bxp.sts$stats) }

  ## Plot boxes
  bxp.data <- bxp(bxp.sts, pars=pars,
                  ylab=ylab,
                  names=names,
                  main=main, cex.main=list(...)$cex.main)
  
  ## add any censored or additional unusual points outside the
  ## boxes
  if(length(bxp.sts$out.endpt) > 0) {
    
    if(logscale) bxp.sts$out.endpt <- log10(bxp.sts$out.endpt)

    with(bxp.sts,
         {
           if(length(out.groupnumber[out.status==1]) > 0) {
             points(out.groupnumber[out.status==1],
                    out.endpt[out.status==1], pch=1)
           }
           if(length(out.groupnumber[out.status==0]) > 0) {
             text(out.groupnumber[out.status==0],
                  out.endpt[out.status==0], labels=">",
                  srt=90)
           }
           if(length(out.groupnumber[out.status==2]) > 0) {
             text(jitter(out.groupnumber[out.status==2]),
                  out.endpt[out.status==2], labels="<",
                  srt=90)
           }
           
         }
         )
  }

  ## also place any censored values that occur at the
  ## box edges (25th and 75th estimated percentiles)
  ## or whiskers
  if(length(bxp.sts$edges.endpt) > 0) {
    
    if(logscale) { bxp.sts$edges.endpt <- log10(bxp.sts$edges.endpt) }

    with(bxp.sts,
         {
           if(length(edges.groupnumber[edges.status==0]) > 0) {
             text(edges.groupnumber[edges.status==0],
                  edges.endpt[edges.status==0], labels=">",
                  srt=90)
           }
           if(length(edges.groupnumber[edges.status==2]) > 0) {
             text(edges.groupnumber[edges.status==2],
                  edges.endpt[edges.status==2], labels="<",
                  srt=90)
           }
         }
         )
  }
  
  invisible()
}


descriptive.censoreddata <- function(dfru, logscale, digits=NULL, offset=NULL,
                                     zeroscore=NULL,  ...) {
  ## PURPOSE: Experimental Function for descriptive summary when censored data
  ## is present. Analogous to boxplotcensoreddata function.
  ##
  ## NOTE: Currently this is designed only for use within the descriptiveTable
  ## method for one factor data. Therefore minimal input checking is done.

  options(warn=-1)
  scfit <-  survfit(survival::Surv(endpt1, endpt2,
                                   status, "interval") ~ -1 + grpf,
                    data=dfru)
  
  summ.scfit <- summary(scfit, censored=TRUE)
  options(warn=0)

  ## Use estimated survival distribution function
  ## values to get quantiles
  scfit.dfr <- with(summ.scfit,
                    ## See comment on scfit.dfr assignment in 
                    ## boxplotcensoreddata() for details.
                    data.frame(grpf=factor(sub ("^grpf=", "", strata),
                                 levels=levels(dfru$grpf)),
                               endpt=time,
                               sdf.prob=surv,
                               cdf.prob=1-surv))
  
  scfit.dfr <- merge(dfru, scfit.dfr, by=c("grpf","endpt"))
  scfit.dfr <- scfit.dfr[order(scfit.dfr$grpf), ]
  grpnames <-  levels(dfru$grpf)

  comps <- lapply(split(scfit.dfr, scfit.dfr$grpf),
                  function(x, logscale=logscale, offset=offset,
                           zeroscore=zeroscore) {
                    y <- x[x$status==1, ]
                    
                    thegrpf <- as.character(unique(x$grpf))
                    ## Follow survfit conventions for estimating
                    ## the median and other quantiles. 
                    ##
                    p50 <- qminmin(x$sdf.prob, x$endpt, 0.50)
                    p25 <- qminmin(x$sdf.prob, x$endpt, 0.25)
                    p75 <- qminmin(x$sdf.prob, x$endpt, 0.75)

                    if(!is.null(offset)) {
                      p25 <- p25 - offset
                      p50 <- p50 - offset
                      p75 <- p75 - offset
                    }

                    ## uncensored and left censored on lower end
                    which.obsmin <- which.min(x$endpt)
                    obsmin.status <- x$status[which.obsmin]
                    if(length(obsmin.status)==1 && obsmin.status==0)  {  ## right censored
                      obsmin <-  NA  ## indeterminate situation
                    }
                    else {
                      obsmin <- min(x$endpt[is.element(x$status, c(1,2))])
                      ## mark only possible left censored observation (status=2)
                      if(max(obsmin.status)==2) {
                        if(!is.null(offset)) { obsmin <- obsmin - offset }
                        if(!is.null(zeroscore)) { obsmin <- 0 }
                        obsmin <- paste("<", obsmin, sep="")
                      }
                    }
                    
                    ## uncensored and right censored on upper end
                    which.obsmax <- which.max(x$endpt)
                    obsmax.status <- x$status[which.obsmax]
                    if(length(obsmax.status)==1 && obsmax.status==2)  {  ## left censored 
                      obsmax <-  NA  ## indeterminate situation
                    }
                    else {
                      obsmax <- max(x$endpt[is.element(x$status, c(1,0))])
                      ## mark only possible right censored observation (status=0)
                      if(min(obsmax.status)==0) {
                        if(!is.null(offset)) { obsmax <- obsmax - offset }
                        obsmax <- paste(">", obsmax, sep="")
                      }
                    }
                    
                    allendpt <- x$endpt
                    allstatus <- x$status

                    allcensendpt <- allendpt[allstatus!=1]
                    allcensstatus <- allstatus[allstatus!=1]
                    
                    thestats <- c(length(allendpt), length(allcensendpt),
                                  length(allendpt) -  length(allcensendpt),
                                  obsmin, p25, p50, p75, obsmax)
                    
                    if(all(allstatus==1)) {
                      thestats["Mean"]==mean(allendpt)
                      if(!is.null(offset)) {
                        thestats[9] <- thestats[9] - offset
                      }
                      thestats["StdDev"]==sd(allendpt)
                      thestats["StdErr"]==sd(allendpt)/sqrt(length(allendpt))
                      
                      if(logscale) {
                        thestats[c("GeoMean","SEGeoMean")] <-
                          c(geoMean(allendpt),
                            geoMean(allendpt)*stndErr(log(allendpt)))
                        if(!is.null(offset)) {
                          thestats[12] <- thestats[12] - offset
                        }
                      }                     
                      
                    }
                    else { ## any censoring at all
                      if(logscale) {
                        thestats[c("Mean","StdDev","StdErr",
                                   "GeoMean","SEGeoMean")] <- rep(NA, 5)
                      }
                      else {
                        thestats[c("Mean","StdDev","StdErr")] <- rep(NA, 3)
                      }
                    }
                    
                    return(thestats)
                    
                  }, logscale=logscale, offset=offset, zeroscore=zeroscore)

  comps <-  as.data.frame(do.call(rbind, comps), stringsAsFactors=FALSE)
  names(comps)[1:11] <- c("n","ncensored","ncomplete","Min",
                          "25%ile","Median",
                          "75%ile","Max",
                          "Mean","StdDev","StdErr")
  if(logscale) {
    names(comps)[12:13] <- c("GeoMean", "SEGeoMean")
  }

  return(comps)
}


## Descriptive Table
setClass("cgOneFactorDescriptiveTable",
         representation(contents="data.frame", settings="list"),
         prototype(contents=data.frame(), settings=list()))

setMethod("descriptiveTable", "cgOneFactorData",
          descriptiveTable.cgOneFactorData <-
          function(data, display="print", ...) {
            ##
            ## PURPOSE: Create a table of quantiles, summary statistics of the data's
            ## individual groups.
            ##
            ## Note that no rounding is done in the calculations
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

            if(!data@has.censored) {
              dlist <-  with(data@dfru, split(endpt, grpf))
              
              comps <- sapply(dlist,
                              function(x) {
                                x <- x[!is.na(x)]
                                c(length(x),min(x),
                                  quantile(x,c(0.25,0.50,0.75)),
                                  max(x),mean(x),sd(x),
                                  stndErr(x) )
                              } )
              dimnames(comps)[[1]] <- c("n","Min","25%ile","Median","75%ile","Max",
                                        "Mean","StdDev","StdErr")
              if(logscale) {
                comps <- rbind(comps, sapply(dlist,
                                             function(x) {
                                               x <- x[!is.na(x)]
                                               c(geoMean(x),
                                                 geoMean(x)*
                                                 stndErr(log(x)),
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
              }

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
              comps <- as.data.frame(t(comps), stringsAsFactors=FALSE)
            }
            
            else {  ## Censored Data
              comps <- descriptive.censoreddata(data@dfru,
                                                logscale=logscale,
                                                settings$digits,
                                                settings$offset,
                                                settings$zeroscore)
            }

            row.names(comps) <- settings$grpnames

            returnObj <- new("cgOneFactorDescriptiveTable",
                             contents=comps,
                             settings=data@settings)
            
            if(display=="print") {
              print(returnObj)
            }
            else if(display=="show") {
              show(returnObj)
            }
            ## else display=="none"
            
            invisible(returnObj)
          }

          )


## Descriptive Table
setMethod("print", "cgOneFactorDescriptiveTable",
          print.cgOneFactorDescriptiveTable <-
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

            ## Accomodate potential of censored representations
            contents <- object@contents
            out <- data.frame(sapply(contents, fround.charcens,
                                     digits),
                              stringsAsFactors=FALSE, check.names=FALSE)
            out[,"n"] <- fround(as.numeric(out[, "n"]), 0)
            if(names(out)[2]=="ncensored") {
              out[, c("ncensored",
                      "ncomplete")] <- sapply(out[,
                                                  c("ncensored","ncomplete")],
                                              fround.charcens, 0)
            }
            row.names(out) <- row.names(contents)
            
            curwidth <- getOption("width")
            on.exit(options(width=curwidth), add=TRUE)
            if(curwidth < 500) { options(width=500)}
            print(out, quote=FALSE)
            invisible()
          })

setMethod("show", "cgOneFactorDescriptiveTable",
          show.cgOneFactorDescriptiveTable <- function(object) showDefault(object))

## Rates Graph for Survival Data, or also, special case
## of Empirical CDF when no censoring is involved

setMethod("kmGraph", "cgOneFactorData",
          kmGraph.cgOneFactorData <-
          function(data, cgtheme=TRUE, distfcn="survival",
                   ylab=NULL, title=NULL,
                   ...) {
            ##
            ## PURPOSE: Graph Nonparametric (default
            ## Kaplan-Meier) survival-type or empirical
            ## cumulative distribution function curves for the groups
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names=c("logscale", "ticklabels"))
            validBoolean(cgtheme)
            
            distfcn <- validArgMatch(distfcn, c("survival","cumulative"))
            if(distfcn=="survival") {
              if(is.null(ylab)) ylab <- "Probability of Survival"
              if(is.null(title)) title <- "Survival Curve Estimates"
            }
            else { ## distfcn=="cumulative"
              if(is.null(ylab)) ylab <- "Cumulative Probability"
              if(is.null(title)) title <- "Cumulative Distribution Function Estimates"
            }
            ## 
            dfru <- data@dfru
            settings <- data@settings
            
            digits <- settings$digits
            offset <- settings$addconstant
            zeroscore <- settings$zeroscore
            
            grpf <- dfru$grpf
            endpt <- dfru$endpt

            if(!is.null(dfru$status)) {
              status <- dfru$status
            }
            else {
              status <- rep(1, length(grpf)) 
            } ## handle all uncensored data case 
            
            numberofgrps <- length(levels(grpf))

            endptlabel <- makeEndptLabel(settings$endptname,
                                         settings$endptunits)

            ticklabelsarg <- getDotsArgName(dots, "ticklabels")
            if(!is.na(ticklabelsarg)) {
              ticklabels <- eval(parse(text=paste("dots$", ticklabelsarg, sep="")))
              validList(ticklabels, names=c("mod","marks"),
                        argname="ticklabels")
            }
            else {
              ticklabels <- NULL
            }
            
            options(warn=-1)
            if(all(dfru$status==1)) { ## no censored values
              thesurvfit <-  survfit(survival::Surv(endpt, status) ~ -1 + grpf,
                                     data=dfru)
            }
            else { ## some censored values
              thesurvfit <-  survfit(survival::Surv(endpt1, endpt2,
                                                    status, "interval") ~ -1 + grpf,
                                     data=dfru)
            }
            options(warn=0)

            ## See comment on scfit.dfr assignment in 
	    ## boxplotcensoreddata() for details.
            stratanames <- sub ("^grpf=", "", names(thesurvfit$strata))
            dfr.survfit <- with(thesurvfit,
                                rbind(data.frame(grp=rep(stratanames,
                                                   times=strata),
                                                 endpt=time,
                                                 prob=surv),
                                      data.frame(grp=stratanames,
                                                 endpt=rep(0,
                                                   length(stratanames)),
                                                 prob=rep(1.0,
                                                   length(stratanames)))))
            ## Prevent survfit's forcing of zero on the axis if there are no zeroes
            if(all(endpt > 0)) {
              dfr.survfit <- dfr.survfit[dfr.survfit$endpt > 0, ]
            }
            
            dfr.survfit$grp <- as.character(dfr.survfit$grp)
            dfr.survfit$grpf <- with(dfr.survfit, factorInSeq(grp))

            dfr.survfit <- with(dfr.survfit,
                                dfr.survfit[order(grpf, endpt), ])

            cgDevice(cgtheme=cgtheme)
            trellispanelstg <- trellis.par.get("clip")$panel
            trellis.par.set("clip", list(panel="off"))
            on.exit(trellis.par.set("clip", list(panel=trellispanelstg)),
                    add=TRUE)
            trellisparstg2 <- trellis.par.get("layout.widths")$ylab
            trellis.par.set("layout.widths", list(ylab=2))
            on.exit(trellis.par.set("layout.widths", list(panel=trellisparstg2)),
                    add=TRUE)
            trellisparstg3 <- trellis.par.get("axis.components")
            trellis.par.set("axis.components",
                            list(left=list(pad1=0.2, pad2=2),
                                 bottom=list(pad1=0.2, pad2=1),
                                 top=list(pad1=0.2, pad2=1),
                                 right=list(pad1=0.2, pad2=2)
                                 ))
            on.exit(trellis.par.set("axis.components", trellisparstg3),
                    add=TRUE)            

            if(!is.null(offset)) {
              endpt <- endpt - offset
            }
            else if(!is.null(zeroscore)) {
              endpt[endpt==min(endpt)] <- 0
            }

            if(distfcn=="cumulative") {
              dfr.survfit$prob <- 1 - dfr.survfit$prob
            }

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

            if(logscale) {
              dfr.survfit$log10endpt <- log10(dfr.survfit$endpt)
              thegraph <- xyplot(prob ~ log10endpt, groups=grpf,
                                 data=dfr.survfit, 
                                 subscripts=TRUE,
                                 digits=digits,
                                 endptlabel=endptlabel,
                                 prepanel=function(x, y, ...) {
                                   prepanel.xYplot(x, y, ...)
                                 },
                                 panel=function(x, y, subscripts, groups, ...)
                                 {
                                   tickmarks <- setupAxisTicks(10^x, axis="x",
                                                               logscale=TRUE,
                                                               grid=TRUE,
                                                               digits=digits,
                                                               offset=offset)
                                   if(!is.null(ticklabels)) {
                                     validList(ticklabels, names=c("mod","marks"),
                                               argname="ticklabels")
                                     tickmarks <- makeTickMarks(ticklabels,
                                                                tickmarks,
                                                                offset=offset)
                                   }
                                   panel.axis(side="bottom",
                                              at=log10(tickmarks),
                                              labels=names(tickmarks),
                                              tck=0.20, text.cex=1,
                                              rot=0,
                                              outside=TRUE)
                                   
                                   panel.grid(h=-1, v=-1)
                                   panel.xYplot(x, y, subscripts, groups,
                                                label.curves=list(method="offset", 
                                                  cex=0.9, adj=0), 
                                                col=cgLineColors[1:numberofgrps],
                                                lty=0, type="l", 
                                                ...)
                                   panel.superpose(x, y, subscripts, groups,
                                                   col=cgLineColors[1:numberofgrps],
                                                   lty=1, lwd=2, type="s", ...)
                                   minmaxTicks(10^x, theaxis="x", logscale=TRUE,
                                               digits=digits,
                                               grid=TRUE,
                                               offset=offset, cex=0.6)
                                 },
                                 xlab=endptlabel,
                                 ylab=ylab,
                                 xlim=rangeExtend(range(dfr.survfit$log10endpt)),
                                 scales=list(
                                   x=list(
                                     alternating=1,
                                     at=range(dfr.survfit$log10endpt),
                                     labels=NULL,
                                     tck=0,
                                     axs="i"),
                                   y=list(alternating=1, tck=0.15)),
                                 main=list(
                                   label=paste(title, "\n",
                                     settings$analysisname, sep=""), cex=1.1)
                                 )
              

            }
            else {
              thegraph <- xyplot(prob ~ endpt, groups=grpf,
                                 data=dfr.survfit, 
                                 subscripts=TRUE,
                                 digits=digits,
                                 endptlabel=endptlabel,
                                 prepanel=function(x, y, ...) {
                                   prepanel.xYplot(x, y, ...)
                                 },
                                 panel=function(x, y, subscripts, groups, ...) {
                                   tickmarks <- setupAxisTicks(endpt, axis="x",
                                                               logscale=FALSE,
                                                               grid=TRUE,
                                                               offset=offset)
                                   if(!is.null(ticklabels)) {
                                     validList(ticklabels, names=c("mod","marks"),
                                               argname="ticklabels")
                                     tickmarks <- makeTickMarks(ticklabels,
                                                                tickmarks,
                                                                offset=offset)
                                   }
                                   panel.axis(side="bottom",
                                              at=tickmarks,
                                              labels=names(tickmarks),
                                              tck=0.20, text.cex=1,
                                              rot=0,
                                              outside=TRUE)
                                   
                                   panel.grid(h=-1, v=-1)
                                   panel.xYplot(x, y, subscripts, groups,
                                                label.curves=list(method="offset", 
                                                  cex=0.9, adj=0), 
                                                col=cgLineColors[1:numberofgrps],
                                                lty=0, type="l", 
                                                ...)
                                   panel.superpose(x, y, subscripts, groups,
                                                   col=cgLineColors[1:numberofgrps],
                                                   lty=1, lwd=2, type="s", ...)
                                   minmaxTicks(endpt, theaxis="x", logscale=FALSE,
                                               digits=digits,
                                               grid=TRUE,
                                               offset=offset, cex=0.6)
                                 },
                                 xlab=endptlabel,
                                 ylab=ylab,
                                 xlim=rangeExtend(range(endpt)),
                                 scales=list(
                                   x=list(
                                     alternating=1,
                                     at=range(endpt),
                                     labels=NULL,
                                     tck=0,
                                     axs="i"),
                                   y=list(alternating=1, tck=0.15)),
                                 main=list(
                                   label=paste(title, "\n",
                                     settings$analysisname, sep=""), cex=1.1)
                                 )
              
            }
            print(thegraph)

            if(settings$stamps) graphStampCG()
            
            invisible()
            
          }
          )



