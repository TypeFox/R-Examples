## $Id: p05SampleSizeOneFactorData.R 6053 2015-02-22 20:23:45Z bpikouni $ 
## One-Factor Unpaired Groups Case

## Sample Size Calculations for One-Factor Unpaired Groups Data

samplesize <- function(sigmaest, endptscale,
                       planningname="", endptname="",
                       ngrps=2, direction="increasing",
                       mmdvec, dendf, ncp=NULL,
                       alpha=0.05, power=0.80,
                       nmax=1000, display="print", ...) {
  ##
  ## PURPOSE: Compute sample sizes based on variability and the other usually
  ## needed specifications. Idea of minimum difference with global
  ## F test is implemented.
  ##
  ## Input argument checking
  validNumeric(ngrps, positive=TRUE, integer=TRUE)
  if(ngrps < 2) {
    stop(cgMessage("The ngrps argument needs to be an integer of 2 or greater."))
  }
  validNumeric(sigmaest, positive=TRUE)
  validCharacter(planningname)
  endptscale <- validArgMatch(endptscale, c("log", "original"))
  direction <- validArgMatch(direction, c("increasing", "decreasing"))
  validNumeric(mmdvec, positive=TRUE)
  validAlpha(alpha)
  validPower(power)
  validDenDf(dendf)
  validNcp(ncp)
  validNumeric(nmax, positive=TRUE, integer=TRUE)
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

  display <- validArgMatch(display, c("print","none","show"))
  dots <- list(...)
  validDotsArgs(dots, names=c(""))

  nsolution <- vector("numeric", length=length(mmdvec))
  numdf <- ngrps - 1
  
  ## main loop to calculate sample size
  for(i in seq(along=mmdvec)) {
    n <- 2 # initialization
    mmdveci <- mmdvec[i]

    if(endptscale=="log") {
      if(direction=="decreasing") {
        sgn <- -1
      }
      else { sgn <- 1 } 
      mmdveci <- log(sgn*mmdveci/100 + 1)
    }

    repeat {
      dendf.result <- do.call("dendf", list(n, ngrps))
      
      if(is.null(ncp)) {
        ## canonical case of one-factor unpaired groups
        ncp.result <- ( n * (mmdveci)^2 ) / ( 2 * sigmaest^2 )
      }
      else {
        ## take the argument for specialized ncp function
        ## which requires this specific order of argument
        ## specification: sigmaest, mmdveci, and n
        ncp.result <- do.call("ncp", list(sigmaest, mmdveci, n))
      }
      pwr <- 1 - pf(qf(1 - alpha, numdf, dendf.result), 
                    numdf, dendf.result, ncp.result)
      if (pwr > power || n >= nmax) {
        nfinal <- n
        break 
      }
      else {
        n <- n + 1
      }
    }
    nsolution[i] <- nfinal 
  }

  if(any(nsolution==nmax)) {
    cat(cgMessage(" With the variability set at ",
                  signif(sigmaest, 4),
                  ", the nmax threshold specified at ",
                  nmax, 
                  "was reached for at least ",
                  "one of the specified differences",
                  "\n\n",
                  warning=TRUE))
  }
  
  computedss <- cbind(mmdvec, n=nsolution, N=ngrps*nsolution)
  
  dimnames(computedss)[[2]] <- c("mmd","n","N")
  n <- computedss[,"n"]
  N <- computedss[,"N"]

  if(display=="print") {
    fmtcomputedss <- as.data.frame(computedss)
    fmtcomputedss$mmd <- fround(fmtcomputedss$mmd, getNumDigits(fmtcomputedss$mmd))
    rownames(fmtcomputedss) <- fmtcomputedss$mmd
    fmtcomputedss$mmd <- NULL
    names(fmtcomputedss) <- c("n per group", "N Total")
    fmtcomputedss$"n per group" <- makeCensored(fmtcomputedss$"n per group",
                                                nmax)
    fmtcomputedss$"N Total" <- makeCensored(fmtcomputedss$"N Total",
                                            ngrps*nmax)

    cat(paste("Sample Size Table",
              if(planningname!="") paste("for", planningname), "\n"))
    if(is.expression(endptname) || (endptname!="")) { 
      cat(paste("Endpoint:", endptname, "\n")) 
    }
    
    ## Taken from base:::chartr help file
    .simpleCap <- function(x) {
      s <- strsplit(x, " ")[[1]]
      paste(toupper(substring(s, 1,1)), substring(s, 2),
            sep="", collapse=" ")
    }
    
    diffmetric <- "Differences"
    if(endptscale=="log") {
      diffmetric <- paste("Percent", .simpleCap(direction), diffmetric)
    }
    else {
      diffmetric <- paste(.simpleCap(direction),  diffmetric)
    }
    cat(diffmetric, "\n")
    
    cat(paste(round(100*power,0),"% Power and ",
              round(100*alpha,0),
              "% Significance Level\n", sep=""))
    
    cat("Variability Estimate ",
        if(endptscale=="log") {
          "(Log scale) "
        },
        "of ",
        signif(sigmaest, 4),
        "\n",
        ngrps, " Groups",
        "\n",
        sep="")

    curwidth <- getOption("width")
    on.exit(options(width=curwidth), add=TRUE)
    if(curwidth < 500) { options(width=500) }

    print(fmtcomputedss, quote=FALSE)
    
  }
  else if(display=="show") {
    showDefault(computedss)
  }
  ## else show nothing
  invisible(computedss)
}

samplesizegraph <- function(sstable, 
                            Nscale="log",
                            mmdscale="log",
                            difftype,
                            direction,
                            analysisname="",
                            endptname="",
                            alpha=0.05,
                            power=0.80,
                            sigmaest,
                            nmax,
                            Nticklabels=NULL,
                            mmdticklabels=NULL,
                            ylab=NULL, ylabright=NULL,
                            titlestamp=TRUE,
                            explanation=TRUE, ...) {
  ##
  ## PURPOSE: Graph Sample Size as a function of Difference
  ##
  ## NOTE: argument sstable needs to be a dataframe of proper format
  ##
  difftype <- validArgMatch(difftype, c("percent","amount","simple"))
  direction <- validArgMatch(direction, c("increasing","decreasing"))
  Nscale <- validArgMatch(Nscale, c("log","original"))
  mmdscale <- validArgMatch(mmdscale, c("log","original"))
  validAlpha(alpha)
  validBoolean(titlestamp)
  validBoolean(explanation)
  dots <- list(...)
  validDotsArgs(dots, names=c(""))
  
  if(length(mmdvec <- sstable[, 1])==1) {
    stop(cgMessage("No graph is produced since only one",
                   "difference value was specified in the sample size",
                   "table"))
  }

  options(warn=-1)
  curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010)
  options(warn=0)
  on.exit(par(curpar))

  n <- sstable[, 2]
  N <- sstable[, 3]
  numberofgrps <- round((N/n)[1], 0)

  N.scaled <- scaleVar(N, endptscale=Nscale)
  N.logscale <- if(Nscale=="log") TRUE else FALSE

  mmd.logscale <- if(mmdscale=="log") TRUE else FALSE
  logandpct <- if(mmd.logscale &&  difftype=="percent") TRUE else FALSE
  mmd.scaled <- scaleVar(mmdvec, endptscale=mmdscale,
                         percent=logandpct)

  parmar <- par(mar=c(5, 4, 4, 3) + 0.1)
  curpar$mar <- parmar$mar

  plot(mmd.scaled, N.scaled, type="n", pch=1, 
       xlab="",
       ylab=if(is.null(ylab)) paste("Total Sample Size from",
         numberofgrps, "groups") else ylab,
       axes=FALSE,
       xlim=rangeExtend(range(mmd.scaled),
         pct=list(minside=6, maxside=0)))
  grid(lty=1)
  lines(mmd.scaled, N.scaled, type="b", pch=1)

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
                    ": ", if(numberofgrps > 2) {"Minimum "} ,"Detectable",
                    " Difference", sep="")

  mtext(side=1, text=xlabchar, line=2)
  if(is.expression(endptname)) {
    mtext(side=1, text=catCharExpr("in ", endptname), line=3)
  }
  else {
    mtext(side=1, text=paste("in", endptname), line=3)
  }
  

  if(titlestamp) {
    title(main=paste("Sample Size Graph\n"), line=2, cex.main=1.1)
  }

  if(explanation) {
    mtext(side=3, line=0.25,
          text=paste(round(100*power,0)," % Power ; ",
            round(100*alpha,0),
            " % Significance Level ; ",
            "Variability Estimate ",
            if(difftype=="percent") {
              "(Log scale) "
            },
            "of ",
            signif(sigmaest, 4), 
            sep=""),
          cex=0.7)
  }
  
  ## Filling in axes tick marks and such.
  names(N.scaled) <- as.character(N)
  N.tickmarks <- setupAxisTicks(N, logscale=N.logscale,
                                digits=0)

  if(!is.null(Nticklabels)) {
    validList(Nticklabels, names=c("mod","marks"),
              argname="Nticklabels")
    N.tickmarks <- makeTickMarks(Nticklabels, N.tickmarks)
  }

  n.ticklabels <- paste("", round(N.tickmarks/numberofgrps, 0))
  ycex <- yTicksCex(N.tickmarks, cexinit=0.8)
  
  axis(2, at=scaleVar(N.tickmarks, endptscale=Nscale),
       labels=names(N.tickmarks), adj=1, cex.axis=ycex, las=1)
  axis(4, at=scaleVar(N.tickmarks, endptscale=Nscale),
       labels=n.ticklabels, adj=0, cex.axis=min(ycex, 0.7), las=1)
  minmaxTicks(N, logscale=N.logscale, digits=0)

  mmdvec.ticks <- if(logandpct) {
    pctToRatio(mmdvec)
  }
  else {
    mmdvec
  }

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
  mtext(side=4, text=if(is.null(ylabright)) paste("Sample size per each of the",
                  numberofgrps, "groups") else ylabright, line=1.8, cex=0.9, adj=0)

  if(any(n >= nmax)) {
    .R. <- TRUE
    coordinates <- Hmisc::largest.empty(mmd.scaled, N.scaled,
                                        width=0.30*(diff(range(
                                          c(min(mmd.scaled),max(mmd.scaled))))),
                                        height=0.25*(diff(range(
                                          c(min(N.scaled),max(N.scaled))))),
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
  
  invisible()
  
}

setClass("cgOneFactorSampleSizeTable",
         representation(ols.sstable="dataframeMatrixOrNULL",
                        rr.sstable="dataframeMatrixOrNULL",
                        settings="list"),
         prototype(ols.sstable=NULL,
                   rr.sstable=NULL, 
                   settings=list()))

setMethod("samplesizeTable", "cgOneFactorFit",
          samplesizeTable.cgOneFactorFit <- 
          function(fit, direction, mmdvec, power=0.80, alpha=0.05,
                   nmax=1000, display="print", ...) {
            ##
            ## PURPOSE:  Compute sample sizes based on variance
            ## estimation from One Factor Fit and other specifications
            ##
            ## Input arguments check
            if(class(fit@uvfit)[1]=="gls" || class(fit@aftfit)[1]=="survreg") {
              stop(cgMessage("There is no samplesizeTable method",
                             "defined for a fitted model that allowed",
                             "unequal variances or censored observations."))
            }
            dots <- list(...)
            validDotsArgs(dots, names=c("model","ngrps"))
            direction <- validArgMatch(direction, c("increasing", "decreasing"))
            validNumeric(mmdvec, positive=TRUE)
            validAlpha(alpha)
            validPower(power)
            validNumeric(nmax, positive=TRUE, integer=TRUE)
            if(nmax < 2) {
              stop(cgMessage("The nmax argument needs to be an integer of 2 or greater."))
            }
            display <- validArgMatch(display, c("print","none","show"))
            modelarg <- getDotsArgName(dots, "model")
            if(!is.na(modelarg)) {
              model <- eval(parse(text=paste("dots$", modelarg, sep="")))
              model <- validArgMatch(model, choices=c("both", "olsonly","rronly"))
            }
            else {
              model <- "both"
            }
            ngrpsarg <- getDotsArgName(dots, "ngrps")
            if(!is.na(ngrpsarg)) {
              ngrps <- eval(parse(text=paste("dots$", ngrpsarg, sep="")))
            }
            else {
              ngrps <- 2
            }

            validNumeric(ngrps, positive=TRUE, integer=TRUE)
            if(ngrps < 2) {
              stop(cgMessage("The ngrps argument needs to be an integer of 2",
                             "or greater."))
            }
           
            ## initializations
            settings <- fit@settings
            ols <- rr <- FALSE
            ols.sstable <- rr.sstable <- NULL
            ols.sigmaest <- rr.sigmaest <- NULL
            ##
            endptscale <- settings$endptscale
            planningname <- settings$analysisname
            endptname <- settings$endptname
            
            rrfit <- fit@rrfit
            olsfit <- fit@olsfit

            if(class(rrfit)[1]=="rlm" && model!="olsonly") {
              rr <- TRUE
            }
            if(class(olsfit)[1]=="lm" && model!="rronly") {
              ols <- TRUE
              if(!rr) model <- "olsonly"
            }

            dendf <- function(n, ngrps) ngrps*(n-1)

            if(ols) {
              ols.sigmaest <- summary(olsfit)$sigma
              ols.sstable <- samplesize(sigmaest=ols.sigmaest,
                                        endptscale=endptscale,
                                        planningname=planningname,
                                        endptname=endptname,
                                        ngrps=ngrps, direction=direction,
                                        mmdvec=mmdvec,
                                        dendf=dendf,
                                        alpha=alpha, power=power,
                                        nmax=nmax,
                                        display="none")
            }
            if(rr) {
              rr.sigmaest <- summary(rrfit)$stddev
              rr.sstable <- samplesize(sigmaest=rr.sigmaest,
                                       endptscale=endptscale,
                                       planningname=planningname,
                                       endptname=endptname,
                                       ngrps=ngrps, direction=direction,
                                       mmdvec=mmdvec,
                                       dendf=dendf,
                                       alpha=alpha, power=power,
                                       nmax=nmax,
                                       display="none")
            }

            ## add to already defined settings
            settings <- c(settings,
                          list(sigmaest=list(ols=ols.sigmaest, rr=rr.sigmaest),
                               planningname=planningname,
                               ngrps=ngrps,
                               direction=direction,
                               alpha=alpha, power=power, nmax=nmax)
                          )
            
            returnObj <- new("cgOneFactorSampleSizeTable",
                             rr.sstable=rr.sstable,
                             ols.sstable=ols.sstable,
                             settings=settings) 

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


setMethod("print", "cgOneFactorSampleSizeTable",
          print.cgOneFactorSampleSizeTable <-
          function(x, title=NULL, endptname=NULL, ...) {
            ##
            ## PURPOSE: Semi-formatted print version of Sample Size Table
            ##
            ## NOTE: Had to use x as an argument because of the system defined
            ## generic. Would have preferred to use object; hence the first
            ## statement below.
            ##
            object <- x
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names="model")

            modelarg <- getDotsArgName(dots, "model")
            if(!is.na(modelarg)) {
              model <- eval(parse(text=paste("dots$", modelarg, sep="")))
              model <- validArgMatch(model, choices=c("both", "olsonly","rronly"))
            }
            else {
              model <- "both"
            }

            ols <- rr <- FALSE  ## Initializations
            if(!is.null(rr.sstable <- object@rr.sstable) && model!="olsonly") {
              rr <- TRUE
              rr.sstable <- as.data.frame(rr.sstable)
            }
            if(!is.null(ols.sstable <- object@ols.sstable) && (model!="rronly")) {
              ols <- TRUE
              ols.sstable <- as.data.frame(ols.sstable)
            }

            settings <- object@settings
            alpha <- settings$alpha
            power <- settings$power
            endptscale <- settings$endptscale
            ngrps <- settings$ngrps
            nmax <- settings$nmax
            
            if(is.null(title)) {
              title <- paste("Sample Size Table from", settings$planningname) 
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

            ## Taken from base:::chartr help file
            .simpleCap <- function(x) {
              s <- strsplit(x, " ")[[1]]
              paste(toupper(substring(s, 1,1)), substring(s, 2),
                    sep="", collapse=" ")
            }
            
            diffmetric <- "Differences"
            if(settings$endptscale=="log") {
              diffmetric <- paste("Percent", .simpleCap(settings$direction), diffmetric)
            }
            else {
              diffmetric <- paste(.simpleCap(settings$direction),  diffmetric)
            }
            cat(diffmetric, "\n")

            cat(paste(round(100*power,0),"% Power and ",
                      round(100*alpha,0),
                      "% Significance Level\n", sep=""))

            fmtdig <- function(x) {
              x$mmd <- fround(x$mmd, getNumDigits(x$mmd))
              rownames(x) <- x$mmd
              x$mmd <- NULL
              names(x) <- c("n per group", "N Total")
              x$"n per group" <- makeCensored(x$"n per group",
                                              nmax)
              x$"N Total" <- makeCensored(x$"N Total",
                                          ngrps*nmax)
              x
            }

            informSettings <- function(sigmaest) {
              cat("Variability Estimate ",
                  if(endptscale=="log") {
                    "(Log scale) "
                  },
                  "of ",
                  signif(sigmaest, 4),
                  "\n",
                  ngrps, " Groups",
                  "\n",
                  sep="")
            }

            curwidth <- getOption("width")
            on.exit(options(width=curwidth), add=TRUE)
            if(curwidth < 500) { options(width=500) }
            
            if(ols) {
              cat("\nClassical Least Squares Based\n")
              informSettings(settings$sigmaest$ols)
              print(fmtdig(ols.sstable), quote=FALSE)
            }
            
            if(rr) {
              cat("\nResistant & Robust Based\n")
              informSettings(settings$sigmaest$rr)
              print(fmtdig(rr.sstable), quote=FALSE)
            }

            invisible()
          })

setMethod("show", "cgOneFactorSampleSizeTable",
          show.cgOneFactorSampleSizeTable <- function(object) showDefault(object))

setMethod("samplesizeGraph", "cgOneFactorSampleSizeTable",
          samplesizeGraph.cgOneFactorSampleSizeTable <-
          function(sstable, Nscale = "log", mmdscale = "log",
                   ## cgtheme=TRUE, device="single",
                   ...) {
            ##
            ## PURPOSE: create a graph of samplesize calculations
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names=c("model",
                                  "mmdticklabels","Nticklabels",
                                  "cgtheme","device"))

            settings <- sstable@settings
            difftype <- if(settings$endptscale=="log") "percent" else "simple"
            alpha <- settings$alpha
            power <- settings$power
            planningname <- settings$planningname
            endptlabel <- makeEndptLabel(settings$endptname, settings$endptunits)   
            ngrps <- settings$ngrps
            direction <- settings$direction
            nmax <- settings$nmax
            stamps <- settings$stamps
            sigmaest <- settings$sigmaest

            rr.sstable <- sstable@rr.sstable
            ols.sstable <- sstable@ols.sstable

            Nticklabelsarg <- getDotsArgName(dots, "Nticklabels")
            if(!is.na(Nticklabelsarg)) {
              Nticklabels <- eval(parse(text=paste("dots$", Nticklabelsarg, sep="")))
              validList(Nticklabels, names=c("mod","marks"),
                        argname="Nticklabels")
            }
            else {
              Nticklabels <- NULL
            }

            mmdticklabelsarg <- getDotsArgName(dots, "mmdticklabels")
            if(!is.na(mmdticklabelsarg)) {
              mmdticklabels <- eval(parse(text=paste("dots$", mmdticklabelsarg, sep="")))
              validList(mmdticklabels, names=c("mod","marks"),
                        argname="mmdticklabels")
            }
            else {
              mmdticklabels <- NULL
            }

            modelarg <- getDotsArgName(dots, "model")
            if(!is.na(modelarg)) {
              model <- eval(parse(text=paste("dots$", modelarg, sep="")))
              model <- validArgMatch(model, choices=c("both", "olsonly","rronly"))
            }
            else {
              dots$model <- "both"
              model <- "both"
            }

            cgthemearg <- getDotsArgName(dots, "cgtheme")
            if(!is.na(cgthemearg)) {
              cgtheme <- eval(parse(text=paste("dots$", cgthemearg, sep="")))
              validBoolean(cgtheme)
            }
            else {
              dots$cgtheme <- TRUE
              cgtheme <- TRUE
            }
            
            devicearg <- getDotsArgName(dots, "device")
            if(!is.na(devicearg)) {
              device <- eval(parse(text=paste("dots$", devicearg, sep="")))
              device <- validArgMatch(device, c("single","multiple", "ask"))
            }
            else {
              dots$device <- "single"
              device <- "single"
            }
            
            ols <- rr <- FALSE  ## Initializations
            if(!is.null(rr.sstable) && model!="olsonly") {
              rr <- TRUE
            }
            if(!is.null(ols.sstable) && model!="rronly") {
              ols <- TRUE
            }


            thetitle <- "Sample Size Graph"

            if(rr && ols && is.element(model, "both") && device=="single") {
              all.dfr <- as.data.frame(rbind(ols.sstable, rr.sstable))
              all.dfr$typef <- factorInSeq(c(rep("Classical",
                                                 nrow(ols.sstable)),
                                             rep("Resistant & Robust",
                                                 nrow(rr.sstable))))
              thetitle <- paste(thetitle, "s", sep="")
              
              cgDevice(cgtheme=cgtheme)
              trellispanelstg <- trellis.par.get("clip")$panel
              trellis.par.set("clip", list(panel="off"))
              on.exit(trellis.par.set("clip", list(panel=trellispanelstg)),
                      add=TRUE)
              trellisparstg2 <- trellis.par.get("layout.widths")$ylab
              trellis.par.set("layout.widths", list(ylab=3))
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
              
              all.N <- all.dfr$N
              all.mmd <- all.dfr$mmd

              N.logscale <- if(Nscale=="log") TRUE else FALSE
              all.N.scaled <- scaleVar(all.N, endptscale=Nscale)

              mmd.logscale <- if(mmdscale=="log") TRUE else FALSE
              logandpct <- if(mmd.logscale &&  difftype=="percent") TRUE else FALSE
              all.mmd.scaled <- scaleVar(all.mmd, endptscale=mmdscale,
                                         percent=logandpct)

              grid.newpage()
              lvp <- viewport(x=0, width=unit(1, "npc") - unit(2, "lines"),
                              name="lvp", just="left")
              tvp <- viewport(x=1, width=unit(2, "lines"),
                              name="tvp", just="right")
              
              thegraph <- xyplot(all.N.scaled  ~ all.mmd.scaled | typef,
                                 data=all.dfr,
                                 panel=function(x, y, ...) {
                                   panel.grid(h=-1, v=-1)
                                   panel.xyplot(x, y, type="b", pch=1,
                                                col="black", ...)
                                   
                                   N.tickmarks <- setupAxisTicks(all.N,
                                                                 logscale=N.logscale,
                                                                 digits=0,
                                                                 grid=TRUE,
                                                                 ycex=0.7)
                                   if(!is.null(Nticklabels)) {
                                     validList(Nticklabels, names=c("mod","marks"),
                                               argname="Nticklabels")
                                     N.tickmarks <- makeTickMarks(Nticklabels, N.tickmarks)
                                   }
                                   
                                   n.ticklabels <- paste("", round(N.tickmarks/ngrps, 0))
                                   ycex <- yTicksCex(N.tickmarks,
                                                     cexinit=0.7,
                                                     cexthreshold = 0.7,
                                                     grid=TRUE)
                                   
				   all.mmd.ticks <- if(logandpct) {
                                     pctToRatio(all.mmd)
                                   }
                                   else {
                                     all.mmd
                                   }

                                   mmdDigits <- getNumDigits(all.mmd)
                                   mmd.tickmarks <-
                                     setupAxisTicks(
                                                    all.mmd.ticks,
                                                    axis="x",
                                                    logscale=mmd.logscale,
                                                    digits=mmdDigits,
                                                    ratio=logandpct,
                                                    percent=logandpct,
                                                    remticks=TRUE,
                                                    grid=TRUE,
                                                    xcex=0.70)

                                   if(length(mmd.tickmarks) < 3) {
                                     mmd.tickmarks <- setupAxisTicks(
                                                                     all.mmd.ticks,
                                                                     axis="x",
                                                                     logscale=mmd.logscale,
                                                                     digits=getNumDigits(all.mmd.ticks),
                                                                     ratio=logandpct,
                                                                     percent=logandpct,
                                                                     remticks=TRUE,
                                                                     grid=TRUE,
                                                                     xcex=0.70)
                                   }

                                   mmd.tickmarks <- mmd.tickmarks[names(mmd.tickmarks)!="0"]
                                   
                                   if(!is.null(mmdticklabels)) {
                                     validList(mmdticklabels, names=c("mod","marks"),
                                               argname="mmdticklabels")
                                     mmd.tickmarks <- makeTickMarks(mmdticklabels, mmd.tickmarks,
                                                                    percent=logandpct)
                                   }
                                   xcex <- xTicksCex(mmd.tickmarks,
                                                     cexinit=0.7,
                                                     cexthreshold = 0.7,
                                                     grid=TRUE)

                                   if(panel.number()==1) {
                                     panel.axis(side="left",
                                                at=scaleVar(N.tickmarks,
                                                  endptscale=Nscale),
                                                labels=names(N.tickmarks),
                                                tck=0.15, text.cex=ycex,
                                                rot=0,
                                                outside=TRUE)
                                   }
                                   else if(panel.number()==2) {
                                     panel.axis(side="right",
                                                at=scaleVar(N.tickmarks,
                                                  endptscale=Nscale),
                                                labels=n.ticklabels,
                                                tck=0.15, text.cex=ycex,
                                                rot=0,
                                                outside=TRUE)
                                   }

                                   panel.axis(side="bottom",
                                              at=scaleVar(mmd.tickmarks,
                                                mmdscale),
                                              labels=names(mmd.tickmarks),
                                              tck=0.15, rot=0,
                                              text.cex=xcex,
                                              outside=TRUE)
                                   N.range <- range(y)
                                   if(N.logscale) {
                                     N.range <- round(10^N.range, 0)
                                   }
                                   minmaxTicks(N.range,
                                               theaxis="y", logscale=N.logscale,
                                               digits=0, 
                                               grid=TRUE)
                                   if(any(round(N.range/ngrps, 0) >= nmax)) {

                                     .R. <- TRUE
                                     coordinates <-  Hmisc::largest.empty(x, y,
                                                                          width=0.30*(diff(range(
                                                                            c(min(x),
                                                                              max(x))))),
                                                                          height=0.25*(diff(range(
                                                                            c(min(y),
                                                                              max(y))))),
                                                                          grid=TRUE,
                                                                          method="area")
                                     panel.text(x=coordinates$rect$x[4], y=coordinates$rect$y[4],
                                                label=paste("For at least one difference\n",
                                                  "the sample size\n",
                                                  "calculations were truncated\n",
                                                  "at",  paste(nmax, "per group."),
                                                  sep=" "),
                                                col="red",
                                                cex=0.6,
                                                adj=c(0, 1))
                                   }
                                 },
                                 layout=c(2,1), aspect=1,
                                 as.table=TRUE,
                                 xlim=rangeExtend(range(all.mmd.scaled),
                                   pct=list(minside=20, maxside=10)),
                                 xlab=paste(if(difftype=="percent") {
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
                                   ": ", if(ngrps > 2) {"Minimum "} ,"Detectable",
                                   " Difference",
                                   sep=""),
                                 ylab=paste("Total Sample Size from",
                                   ngrps, "groups"),
                                 scales=list(y=list(labels=NULL, tck=0),
                                   x=list(labels=NULL, tck=0)),
                                 main=list(label=paste(thetitle, "\n",
                                             planningname,
                                             sep=""), cex=1.1),
                                 par.strip.text=list(cex=0.7)
                                 )                                 
              
              ##on.exit(lattice:::lattice.setStatus(print.more = FALSE),
              ##        add=TRUE)
              pushViewport(lvp)
              print(thegraph, newpage=FALSE)
              upViewport()
              pushViewport(tvp)
              grid.text(paste("Sample size per each of the",
                              ngrps, "groups"),
                        x = unit(1, "lines"), rot=90, gp=gpar(cex=0.75),
                        just="center")
              
              seekViewport(trellis.vpname("xlab"))
              grid.text(catCharExpr("in ", endptlabel), y = unit(-1, "lines"))
              upViewport(0)
              
              seekViewport(trellis.vpname("main"))
              grid.text(paste(round(100*power,0)," % Power ; ",
                              round(100*alpha,0),
                              " % Significance Level ; ",
                              "Variability Estimate ",
                              if(difftype=="percent") {
                                "(Log scale) "
                              },
                              "of\n",
                              signif(sigmaest$ols, 4),
                              " for Classical and ",
                              signif(sigmaest$rr, 4),
                              " for Resistant & Robust",
                              sep=""),
                        y=unit(-1, "lines"), rot=0, gp=gpar(cex=0.60),
                        just="center")
              upViewport(0)
              usedgrid <- TRUE
            }
            
            else if(((model=="olsonly" && ols) || (ols && !rr && model=="both")) &&
                    device=="single") {
              samplesizegraph(ols.sstable,
                              mmdscale,
                              Nscale,
                              difftype,
                              direction,
                              planningname,
                              endptlabel,
                              alpha,
                              power,
                              sigmaest=sigmaest$ols,
                              nmax=nmax,
                              Nticklabels,
                              mmdticklabels,
                              titlestamp=FALSE,
                              explanation=FALSE)
              title(line=2,
                    main=paste("Sample Size Graph\n",
                      planningname, sep=""), cex=1.1)
              mtext(side=3, line=0.5,
                    text=paste(round(100*power,0)," % Power ; ",
                      round(100*alpha,0),
                      " % Significance Level ; ",
                      "Classical Variability Estimate ",
                      if(difftype=="percent") {
                        "(Log scale) "
                      },
                      "of ",
                      signif(sigmaest$ols, 4),
                      sep=""), cex=0.60)
              usedgrid <- FALSE
              
            }

            else if(((model=="rronly" && rr) || (!ols && rr && model=="both")) &&
                    device=="single") {
              samplesizegraph(rr.sstable,
                              mmdscale,
                              Nscale,
                              difftype,
                              direction,
                              planningname,
                              endptlabel,
                              alpha,
                              power,
                              sigmaest=sigmaest$rr,
                              nmax=nmax,
                              Nticklabels,
                              mmdticklabels,
                              titlestamp=FALSE,
                              explanation=FALSE)
              title(line=2,
                    main=paste("Sample Size Graph\n",
                      planningname, sep=""), cex=1.1)
              mtext(side=3, line=0.5,
                    text=paste(round(100*power,0)," % Power ; ",
                      round(100*alpha,0),
                      " % Significance Level ; ",
                      "Resistant & Robust Variability Estimate ",
                      if(difftype=="percent") {
                        "(Log scale) "
                      },
                      "of ",
                      signif(sigmaest$rr, 4),
                      sep=""), cex=0.60)
              usedgrid <- FALSE
              
            }
            else if(rr && ols &&
                    is.element(device, c("ask","multiple")) &&
                    is.element(model, "both")) { 
              
              device <- validArgMatch(device, c("multiple", "ask"))
              if(device=="ask") {
                op <- par(ask = TRUE)
                on.exit(par(op), add=TRUE)
              }

              samplesizegraph(ols.sstable,
                              mmdscale,
                              Nscale,
                              difftype,
                              direction,
                              planningname,
                              endptlabel,
                              alpha,
                              power,
                              sigmaest=sigmaest$ols,
                              nmax=nmax,
                              Nticklabels,
                              mmdticklabels,
                              titlestamp=FALSE,
                              explanation=FALSE)
              title(line=2,
                    main=paste("Sample Size Graph\n",
                      planningname, sep=""), cex=1.1)
              mtext(side=3, line=0.5,
                    text=paste(round(100*power,0)," % Power ; ",
                      round(100*alpha,0),
                      " % Significance Level ; ",
                      "Classical Variability Estimate ",
                      if(difftype=="percent") {
                        "(Log scale) "
                      },
                      "of ",
                      signif(sigmaest$ols, 4),
                      sep=""), cex=0.60)
              if(stamps) graphStampCG(grid=FALSE)

              if(device=="multiple") {
                ## since we only want dots arguments for trellis.device in
                ## next call we need some housekeeping 
                if(!is.null(dots$model)) dots$model <- NULL
                if(!is.null(dots$device) &&
                   is.character(validArgMatch(dots$device, c("single","multiple","ask")))) {
                  dots$device<- NULL
                }
                if(!is.null(dots$cgtheme)) dots$cgtheme <- NULL
                if(!is.null(dots$mmdticklabels)) dots$mmdticklabels <- NULL
                if(!is.null(dots$Nticklabels)) dots$Nticklabels <- NULL
                
                do.call("cgDevice", c(list(new=TRUE), dots))
                cat(cgMessage("A new graphics device has been generated",
                              "to hold",
                              "SampleSize graph version based on",
                              "the Resistant & Robust estimate.",
                              "The version based on the Classical Least Squares",
                              "estimate is on the previous",
                              "device.\n",
                              warning=TRUE))
              }
              samplesizegraph(rr.sstable,
                              mmdscale,
                              Nscale,
                              difftype,
                              direction,
                              planningname,
                              endptlabel,
                              alpha,
                              power,
                              sigmaest=sigmaest$rr,
                              nmax=nmax,
                              Nticklabels,
                              mmdticklabels,
                              titlestamp=FALSE,
                              explanation=FALSE)
              title(line=2,
                    main=paste("Sample Size Graph\n",
                      planningname, sep=""), cex=1.1)
              mtext(side=3, line=0.5,
                    text=paste(round(100*power,0)," % Power ; ",
                      round(100*alpha,0),
                      " % Significance Level ; ",
                      "Resistant & Robust Variability Estimate ",
                      if(difftype=="percent") {
                        "(Log scale) "
                      },
                      "of ",
                      signif(sigmaest$rr, 4),
                      sep=""), cex=0.60)
              usedgrid <- FALSE
            }
            else {
              stop(cgMessage("The chosen device and model arguments",
                             "are not compatible either with each other",
                             "or with the fitted model(s) in the SampleSize table object.",
                             seeHelpFile("SampleSizeTable")))
            }

            if(stamps) graphStampCG(grid=usedgrid)
            invisible()
            
          }
          )


validDenDf <- function(dendf) {
  if(!is.function(dendf)) {
    stop(cgMessage("The object needs to be a function",
                   "with arguments n and ngrps."))
  }
  if(!identical(sort(names(formals(dendf))),
                c("n","ngrps"))) {
    stop(cgMessage("The dendf function needs to include",
                   "the exact argument names of n and ngrps."))
  }
  return(TRUE)
}

validNcp <- function(ncp) {
  if(is.null(ncp)) return(TRUE)
  if(!is.function(ncp)) {
    stop(cgMessage("The object needs to be a function."))
  }
  theargnames <- names(formals(ncp))
  if(sum(theargnames %in% c("sigmaest", "mmdvec", "n")) < 3) {
    stop(cgMessage("The ncp function needs to include",
                   "the argument names of sigmaest, mmdvec, n."))
  }
  return(TRUE)
}




