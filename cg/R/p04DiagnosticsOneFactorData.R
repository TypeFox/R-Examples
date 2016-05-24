## $Id: p04DiagnosticsOneFactorData.R 4891 2014-02-21 17:49:45Z bpikouni $
## One-Factor Unpaired Groups Case

## Diagnostics methods for One-Factor Unpaired Groups Data

variancegraph <- function(resids, fitteds, 
                          status=NULL,
                          grp=NULL, 
                          desc="", analysisname="",
                          endptname="", trend=NULL,
                          titlestamp=TRUE,
                          trendstamp=FALSE,  ...) {
  ##
  ## PURPOSE: Create a Tukey-inspired plot of the residuals to assess
  ## variances assumption
  ##
  ## NOTE: Incomplete checking is done if censored residuals are
  ## specified. The survival::Surv() convention is presumed as the sole
  ## required data format. 
  ##
  ## Input Argument Checking
  validAtomicVec(resids)
  validAtomicVec(fitteds)
  if(!is.null(status)) {
    validAtomicVec(status)
    if(any(!is.element(status, c(0,1,2)))) {
      stop(cgMessage("The status argument values",
                     "to indicate censoring can only have valid values",
                     "of 0, 1, or 2, to indicate right, no, or left censoring",
                     "consistent with the Surv() function from the survival",
                     "package. (The variancegraph() function cannot",
                     "handle interval censored data.)"))
    }
    if((length(resids) != length(fitteds)) || (length(resids) != length(status))) {
      stop(cgMessage("The resids, fitteds, and status argument values",
                     "must be have equal length."))
    }
    has.censored <- TRUE    
  }
  else {
    if((length(resids) != length(fitteds)) ) {
      stop(cgMessage("The resids and fitteds argument values",
                     "must be have equal length."))
    }
    has.censored <- FALSE
  }
  validBoolean(titlestamp)
  validBoolean(trendstamp)

  ## End input argument handling
  
  options(warn=-1)
  curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010)
  options(warn=0)
  on.exit(par(curpar), add=TRUE)

  ylabchar <- if((is.character(endptname) && endptname!="" ||
                  is.expression(endptname)))  {
    catCharExpr("Absolute Residual in  ",
                endptname)
  }
  else {
    "Absolute Residual"
  }

  tukeyresids <- sqrt(abs(resids))
  if(has.censored) {
    tukeyresids2 <- sqrt(abs(resids))
  }
  else {
    tukeyresids2 <- NULL
  }
  
  if(!is.null(grp)) {
    grpf <- if(is.factor(grp)) grp else as.numeric(factorInSeq(grp))
    if(has.censored) {
      n.uncensored <- sapply(split(status, grpf),
                             function(x) {
                               sum(x[x==1])
                             })
    }

    ## next assignment creates a data frame with smooths
    dfr <- residualgrptrend.helper(tukeyresids, fitteds, tukeyresids2,
                                   status, grpf, desc)
    grpnames <- unique(dfr$grpc)
    numberofgrps <- length(grpnames)

    if((is.null(trend))) {
      trend <- TRUE ## initialization

      if(has.censored) {
        if(any(n.uncensored < 5)) {
          trend <- FALSE
          trendstamp <- FALSE
          warning(cgMessage("No trend line is displayed",
                            "since at least one of the groups",
                            "has less than 5",
                            "complete observations.",
                            "This can be overidden by setting",
                            "trend=TRUE in the call.",
                            warning=TRUE))
        }
      }
    }

    with(dfr,
         {
           jittergrpn <- jitter(grpn)
           plot(jittergrpn, tukeyresids,
                ylab=ylabchar,
                xlab="",
                xlim=c(0, numberofgrps) + 0.5, axes=FALSE,
                type="n")
           if(has.censored) {
             if(sum(status==0) > 0) {
               text((jittergrpn)[status==0],
                    tukeyresids[status==0], labels=">", srt=90,
                    col="darkgray")
             }
             points((jittergrpn)[status==1], tukeyresids[status==1])
             if(sum(status==2) > 0) {
               text((jittergrpn)[status==2],
                    tukeyresids[status==2], labels="<", srt=90,
                    col="darkgray")
               
             }
           }
           else {
             points(jittergrpn, tukeyresids)
           }
           if(trend) {
             lines(grpn, trendtukeyresids)
             
           }
         })

    if(trendstamp && trend) {
      if(!has.censored) {
        mtext(text=paste("smoothed line with x-axis ordered by fitted group means"),
              side=3,
              cex=0.6, col="gray", line=0.2, adj=0.5)
      }
      else {
        mtext(text=paste("line connects estimated group means",
                "with x-axis ordered by fitted group means"),
              side=3,
              cex=0.6, col="gray", line=0.2, adj=0.5)
      }
    }
    
    if(trend & has.censored && any(n.uncensored < 5)) {
      mtext(side=1, line=4,
            text=paste("NOTE: At least one group has less than 5",
              "complete observations, so trends around it may not make",
              "sense."), col="red", cex=0.6)
    }
    
    ## Axes Customization
    grpnameticksettings <- setupGrpNameTicks(grpnames, 1:numberofgrps)
    plotGrpNameTicks((grpnames), settings=grpnameticksettings)

    ## faint gray horizontal dividers between groups
    abline(v=seq(1.5, numberofgrps-0.5, 1), col="gray")
  }
  else {
    plot(fitteds, tukeyresids,
         ylab=ylabchar,
         xlab="",
         axes=FALSE, type="n")
    points(fitteds, tukeyresids)
    if(trend) {
      lines(lowess(fitteds, tukeyresids))
    }
    x.tickmarks <- setupAxisTicks(fitteds, logscale=FALSE)
    axis(1, at=x.tickmarks, labels=names(x.tickmarks), cex.axis=0.8)
  }

  tickmarks <- setupAxisTicks(tukeyresids^2, logscale=FALSE)
  axis(2, at=sqrt(tickmarks), labels=names(tickmarks),
       cex.axis=0.8, adj=1, las=1)
  mtext("square-root spaced", side=2, line=2.25, cex=0.7)
  
  text(x=rep(par("usr")[1], 2), y=range(tukeyresids),
       labels=paste("", signif(range(tukeyresids^2), digits=3)),
       col="blue", adj=0, cex=0.7)
  box()
  
  if(titlestamp) {
    title(main=paste("Variance Graph: ", desc,
            "\n", analysisname, sep=""), line=2, cex.main=1.1)
  }
  
  invisible()
  
}

setMethod("varianceGraph", "cgOneFactorFit",
          varianceGraph.cgOneFactorFit <-
          function(fit,  trend=NULL, cgtheme=TRUE, device="single", ...) {
            ##
            ## PURPOSE: Create a Tukey-inspired plot of the Residuals to assess
            ## variances assumption of homoscedasticity, perhaps after weighting.
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names=c("model"))
            validBoolean(cgtheme)

            rrfit <- fit@rrfit
            olsfit <- fit@olsfit
            aftfit <- fit@aftfit
            uvfit <- fit@uvfit

            modelarg <- getDotsArgName(dots, "model")
            if(!is.na(modelarg)) {
              model <- eval(parse(text=paste("dots$", modelarg, sep="")))
              model <- validArgMatch(model,c("olsonly","rronly","rrwtdonly",
                                             "rrunwtdonly","extended","both"))
            }
            else {
              model <- "both"
            }

            device <- validArgMatch(device, c("single","multiple", "ask"))

            aft <- ols <- rr <- uv <- FALSE  ## Initializations

            if(class(aftfit)[1]=="survreg") {
              aft <- TRUE
              validArgModel(...)
            }
            else if(class(uvfit)[1]=="gls") {
              uv <- TRUE
              validArgModel(...)
            }
            
            if(class(rrfit)[1]=="rlm" && model!="olsonly" && !aft && !uv) {
              rr <- TRUE
            }
            if(class(olsfit)[1]=="lm" && model!="rronly" && !aft && !uv) {
              ols <- TRUE
              if(!rr) model <- "olsonly"
            }

            settings <- fit@settings
            analysisname <- settings$analysisname
            endptlabel <- makeEndptLabel(settings$endptname, settings$endptunits)

            numberofgrps <- length(settings$grpnames)
            thetitle <- "Variance Graph"
            thesmoothmsg <- paste("smoothed line with x-axis ordered by fitted group means")


            if(rr && ols && is.element(model, c("both", "extended"))
               && device=="single") {
              ols.tukeyresids <- sqrt(abs(resid(olsfit)))
              ols.dfr <- residualgrptrend.helper(ols.tukeyresids,
                                                 fitted(olsfit),
                                                 grpf=olsfit$dfru$grpf,
                                                 desc="Classical")
              rrwtd.tukeyresids <- sqrt(abs(sqrt(rrfit$w)*resid(rrfit)))
              rrwtd.dfr <- residualgrptrend.helper(rrwtd.tukeyresids,
                                                   fitted(rrfit),
                                                   grpf=rrfit$dfru$grpf,
                                                   desc="RR Weighted")
              if(model=="both") {
                all.dfr <- rbind(ols.dfr, rrwtd.dfr)
                all.tukeyresids <- c(ols.tukeyresids, rrwtd.tukeyresids)
                thelayout <- c(2, 1)
              }
              else if(model=="extended") {
                rr.tukeyresids <- sqrt(abs(resid(rrfit)))
                rr.dfr <- residualgrptrend.helper(rr.tukeyresids,
                                                  fitted(rrfit),
                                                  grpf=rrfit$dfru$grpf,
                                                  desc="RR Unweighted")
                all.dfr <- rbind(ols.dfr, rrwtd.dfr, rr.dfr)
                all.tukeyresids <- c(ols.tukeyresids, rrwtd.tukeyresids,
                                     rr.tukeyresids)
                thelayout <- c(3, 1)
              }

              thetitle <- paste(thetitle, "s", sep="")
              thesmoothmsg <- gsub("line", "lines", thesmoothmsg)
              
              all.dfr$typef <- factorInSeq(all.dfr$type)
              all.dfr$grpc <- as.character(all.dfr$grpc)
              grplabels <- with(all.dfr, lapply(split(grpc, typef), unique))

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
              
              thegraph <- xyplot(Cbind(tukeyresids, trendtukeyresids)
                                 ~ grpn | typef,
                                 data=all.dfr,
                                 panel=function(x, y, subscripts=TRUE, ...) {
                                   panel.number <- panel.number()
                                   curgrplabels <- grplabels[[panel.number]]

                                   grpnameticksettings <- setupGrpNameTicks(curgrplabels,
                                                                            grplocation=
                                                                            1:length(curgrplabels),
                                                                            cexinit=0.8,
                                                                            cexthreshold=0.5,
                                                                            grid=TRUE)
                                   plotGrpNameTicks(curgrplabels, grpnameticksettings, grid=TRUE)
                                   if(is.element(panel.number,
                                                 c(1, nlevels(all.dfr$typef)))) {
                                     panelside <- ifelse(panel.number==1,
                                                         "left", "right") 
                                     tickmarks <- setupAxisTicks(all.dfr$tukeyresids^2,
                                                                 logscale=FALSE, grid=TRUE)
                                    panel.axis(side=panelside,
                                                at=sqrt(tickmarks),
                                                labels=names(tickmarks),
                                                tck=0.20, text.cex=0.6,
                                                rot=0,
                                                outside=TRUE)
                                   }
                                   if(is.null(trend) || trend) {
                                     panel.lines(x, unwind(attr(y, "other")),
                                                 ...)
                                   }
                                   panel.points(jitter(x), y, ...)
                                   panel.abline(v=seq(1.5, max(x) - 0.5, by=1),
                                                col="lightgray")
                                   panel.text(x=rep(0, 2), y=range(y),
                                              labels=paste("",
                                                c(signif(min(y^2), 3),
                                                  signif(max(y^2), 3))),
                                              col="blue", adj=0, cex=0.6)
                                 },
                                 col="black", xlim=c(0, numberofgrps + 0.5),      
                                 xlab="",
                                 ylab=list(cex=0.7,
                                   label="  "),
                                 layout=thelayout, aspect=1,
                                 scales=list(x=list(at=1:numberofgrps,
                                               labels=NULL, tck=c(0.15, 0), axs="i"),
                                   y=list(labels=NULL, tck=0)),
                                 main=list(label=paste(thetitle, "\n",
                                             fit@settings$analysisname, sep=""), cex=1.1),
                                 par.strip.text=list(cex=0.7))

              print(thegraph, position=c(0,0,0.95,1))

              seekViewport(trellis.vpname("ylab"))
              ylabchar <- catCharExpr("Absolute Residual in Fit of  ", endptlabel)
              grid.text(ylabchar,
                        x = unit(0, "lines"), rot=90, gp=gpar(cex=0.9))
              grid.text("square-root spaced",
                        x = unit(1, "lines"), rot=90, gp=gpar(cex=0.7))

              grid::seekViewport(trellis.vpname("main"))
              grid.text(thesmoothmsg, gp=gpar(col="gray", cex=0.6),
                        y=unit(-1, "lines"))
              upViewport(0)
              
              if(settings$stamps) graphStampCG()

              invisible()
            }
            
            else if(model=="olsonly" && ols && device=="single") {
              variancegraph(resid(olsfit), fitted(olsfit),
                            grp = olsfit$dfru$grpf, desc = "Classical",
                            analysisname=analysisname, 
                            endptname=catCharExpr(" Fit of ", endptlabel),
                            trend=trend,
                            titlestamp = TRUE, trendstamp=TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)

            }

            else if((model=="rronly" || model=="rrwtdonly") && rr
                    && device=="single") {
              variancegraph(sqrt(rrfit$w)*resid(rrfit), fitted(rrfit),
                            grp = olsfit$dfru$grpf,
                            desc = "Resistant & Robust Weighted",
                            analysisname=analysisname, 
                            endptname=catCharExpr(" Fit of ", endptlabel),
                            trend=trend,
                            titlestamp = TRUE, trendstamp=TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)
            }

            else if(model=="rrunwtdonly" && rr && device=="single") {
              variancegraph(resid(rrfit), fitted(rrfit),
                            grp = olsfit$dfru$grpf,
                            desc = "Resistant & Robust Unweighted",
                            analysisname=analysisname, 
                            endptname=catCharExpr(" Fit of ", endptlabel),
                            trend=trend,
                            titlestamp = TRUE, trendstamp=TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)
            }

            else if(rr && ols &&
                    is.element(device, c("ask","multiple")) &&
                    is.element(model, c("both", "extended"))) {
              if(device=="ask") {
                op <- par(ask = TRUE)
                on.exit(par(op), add=TRUE)
              }
              else { ## device=="multiple"
                cat(cgMessage("New graphic devices will be generated to",
                              "hold multiple graphs requested...\n",
                              warning=TRUE))
              }
              variancegraph(resid(olsfit), fitted(olsfit),
                            grp = olsfit$dfru$grpf, desc = "Classical",
                            analysisname=analysisname, 
                            endptname=catCharExpr(" Fit of ", endptlabel),
                            trend=trend,
                            titlestamp = TRUE, trendstamp=TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)

              if(device=="multiple") {
                do.call("cgDevice", c(list(new=TRUE), dots))
              }
              variancegraph(sqrt(rrfit$w)*resid(rrfit), fitted(rrfit),
                            grp = olsfit$dfru$grpf,
                            desc = "Resistant & Robust Weighted",
                            analysisname=analysisname, 
                            endptname=catCharExpr(" Fit of ", endptlabel),
                            trend=trend,
                            titlestamp = TRUE, trendstamp=TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)

              if(model=="extended") {
                if(device=="multiple") {
                  do.call("cgDevice", c(list(new=TRUE), dots))
                }
                variancegraph(resid(rrfit), fitted(rrfit),
                              grp = olsfit$dfru$grpf,
                              desc = "Resistant & Robust Unweighted",
                              analysisname=analysisname, 
                              endptname=catCharExpr(" Fit of ", endptlabel),
                              trend=trend,
                              titlestamp = TRUE, trendstamp=TRUE)
                if(settings$stamps) graphStampCG(grid=FALSE)
              }
            }
            
            else if(uv && device=="single") {
              variancegraph(residuals(uvfit, type="pearson"), fitted(uvfit),
                            grp = uvfit$dfru$grpf,
                            desc = "Unequal Variances",
                            analysisname = analysisname, 
                            endptname=catCharExpr("  , Standardized  ",
                              catCharExpr(" Fit of ", endptlabel), rev=TRUE),
                            trend=trend,
                            titlestamp=TRUE, trendstamp=TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)
            }

            else if(aft && device=="single") {
              variancegraph(residuals(aftfit, type="response"),
                            predict(aftfit),
                            status=aftfit$dfru$status,
                            grp = aftfit$dfru$grpf,
                            desc = "AFT Censored",
                            analysisname=analysisname, 
                            endptname=catCharExpr(" Fit of ", endptlabel),
                            trend=trend,
                            titlestamp = TRUE, trendstamp=TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)
            }

            else {
              stop(cgMessage("The chosen device and model arguments",
                             "are not compatible either with each other",
                             "or with the fitted model(s) in the fit object.",
                             seeHelpFile("varianceGraph")))
            }

            
            
            invisible()
          }
          )            


qqgraph <- function(resids, line=NULL, status=NULL, 
                    desc="", analysisname="",
                    endptname="",  titlestamp=TRUE,
                    ...) {
  ##
  ## PURPOSE: Create a Q-Q Gaussian Plot(s) of the Residuals
  ## Really just a wrapper around stats::qqplot.default
  ## Input Argument Checking
  validAtomicVec(resids)
  validBoolean(titlestamp)
  if(!is.null(status)) {
    validAtomicVec(status)
    if(any(!is.element(status, c(0,1,2)))) {
      stop(cgMessage("The status argument values",
                     "to indicate censoring can only have valid values",
                     "of 0, 1, or 2, to indicate right, no, or left censoring",
                     "consistent with the Surv() function from the survival",
                     "package. (The variancegraph() function cannot",
                     "handle interval censored data.)"))
    }
    if(length(resids) != length(status)) {
      stop(cgMessage("The resids, fitteds, and status argument values",
                     "must be have equal length."))
    }
    has.censored <- TRUE    
  }
  else {
    has.censored <- FALSE
  }
  ## End input argument handling
  
  options(warn=-1)
  curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010)
  options(warn=0)
  on.exit(par(curpar), add=TRUE)

  ylabchar <- if((is.character(endptname) && endptname!="" ||
                  is.expression(endptname)))  {
    catCharExpr("Residual in  ",
                endptname)
  }
  else {
    "Residual"
  }

  xy <- qqnorm(y=resids, main = "",
               xlab = "Standard Gaussian (Normal) Quantile", 
               ylab=ylabchar,
               type="n", axes=FALSE, ...)
  grid(lty=1)
  if(has.censored) {
    x.cens <- xy$x
    y.cens <- xy$y
    if(sum(status==0) > 0) {
      text(x.cens[status==0],
           y.cens[status==0], labels=">", srt=90,
           col="darkgray")
    }

    points(x.cens[status==1], y.cens[status==1])

    if(sum(status==2) > 0) {
      text(x.cens[status==2],
           y.cens[status==2], labels="<", srt=90,
           col="darkgray")
    }
  }
  else {
    points(xy$x, xy$y)
  }


  if(is.null(line)) {
    if(!has.censored) {
      qqline(as.vector(resids)) ## avoid curious residuals.lm() AsIs bug if needed
    }
    ## otherwise do not graph line if there are any censored data residuals,
    ## i.e. line=TRUE needs to be specified
  }
  else if(line) {
    if(!has.censored) {
      qqline(as.vector(resids)) ## avoid curious residuals.lm() AsIs bug if needed
    }
    else {  ## try to see if 25th and 75th quantiles can be estimated
      ## in the presence of censoring
      ## Follow survfit conventions for estimation
      ##
      suppressWarnings( {
        scfit <-  survfit(survival::Surv(resids, resids,
                                         status, "interval") ~ 1)
        summ.scfit <- summary(scfit, censored=TRUE)
        scfit.dfr <- with(summ.scfit,
                          data.frame(resids=time,
                                     sdf.prob=surv))
        
        scfit.dfr <- merge(data.frame(resids=resids), scfit.dfr,
                           by=c("resids"))
        p25 <- with(scfit.dfr, qminmin(sdf.prob, resids, 0.25))
        p75 <- with(scfit.dfr, qminmin(sdf.prob, resids, 0.75))

        if(is.na(p25) || is.na(p75)) {
          warning(cgMessage("No line is displayed",
                            "since at least one of the 25th or",
                            "the 75th percentiles could not",
                            "estimated from the censored data residuals",
                            warning=TRUE))
        }
        else {
          ## following stats::qqline definition
          y.res <- c(p25, p75)
          x.res <- qnorm(c(0.25, 0.75))
          slope.res <- diff(y.res)/diff(x.res)
          int <- y.res[1L] - slope.res * x.res[1L]
          abline(int, slope.res)
        }
        
      })
    }
  }
  tickmarks <- setupAxisTicks(resids, logscale=FALSE)
  axis(2, at=tickmarks, labels=names(tickmarks),
       cex.axis=0.8, adj=1, las=1)
  axis(1, cex.axis=0.8)
  box()

  text(x=rep(par("usr")[1], 2), y=range(resids),
       labels=paste("", signif(range(resids), digits=3)),
       col="blue", adj=0, cex=0.7)

  if(titlestamp) {
    title(main=paste("Quantile-Quantile (QQ) Graph on\nGaussian (Normal)",
            " Distribution: ", desc,
            "\n", analysisname, sep=""), line=1, cex.main=1.1)
  }
  
  invisible()
  
}

setMethod("qqGraph", "cgOneFactorFit",
          qqGraph.cgOneFactorFit <- 
          function(fit, line=NULL, cgtheme=TRUE, device="single", ...) {
            ##
            ## PURPOSE: Create a Q-Q Gaussian Plot(s) of the Residuals
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names=c("model"))
            validBoolean(cgtheme)

            rrfit <- fit@rrfit
            olsfit <- fit@olsfit
            aftfit <- fit@aftfit
            uvfit <- fit@uvfit

            modelarg <- getDotsArgName(dots, "model")
            if(!is.na(modelarg)) {
              model <- eval(parse(text=paste("dots$", modelarg, sep="")))
              model <- validArgMatch(model,c("olsonly","rronly","rrwtdonly",
                                             "rrunwtdonly","extended","both"))
            }
            else {
              model <- "both"
            }
            
            device <- validArgMatch(device, c("single","multiple", "ask"))
            
            rrfit <- fit@rrfit
            olsfit <- fit@olsfit
            aftfit <- fit@aftfit
            uvfit <- fit@uvfit
            
            aft <- ols <- rr <- uv <- FALSE  ## Initializations
            
            if(class(aftfit)[1]=="survreg") {
              aft <- TRUE
              validArgModel(...)
            }
            else if(class(uvfit)[1]=="gls") {
              uv <- TRUE
              validArgModel(...)
            }
            
            if(class(rrfit)[1]=="rlm" && model!="olsonly" && !aft && !uv) {
              rr <- TRUE
            }
            if(class(olsfit)[1]=="lm" && model!="rronly" && !aft && !uv) {
              ols <- TRUE
              if(!rr) model <- "olsonly"
            }

            settings <- fit@settings
            analysisname <- settings$analysisname
            endptlabel <- makeEndptLabel(settings$endptname, settings$endptunits)

            numberofgrps <- length(settings$grpnames)
            thetitle <- "Quantile-Quantile (QQ) Graph on\nGaussian (Normal) Distribution"
            
            residual.helper <- function(resid, type) {
              return(data.frame(resid=resid,
                                type=rep(type, length(resid))))
            }

            if(rr && ols && is.element(model, c("both", "extended"))
               && device=="single") {
              ols.dfr <-  residual.helper(residuals(olsfit), "Classical")
              rrwtd.dfr <- residual.helper(sqrt(rrfit$w)*resid(rrfit),
                                           "RR Weighted")
              
              if(model=="both") {
                all.dfr <- rbind(ols.dfr, rrwtd.dfr)
                all.resid <- c(ols.dfr$resid, rrwtd.dfr$resid)
                thelayout <- c(2, 1)
              }
              else if(model=="extended") {
                rrunwtd.dfr <- residual.helper(resid(rrfit),
                                               "RR Unweighted")
                all.dfr <- rbind(ols.dfr, rrwtd.dfr, rrunwtd.dfr)
                all.resid <- c(ols.dfr$resid, rrwtd.dfr$resid,
                               rrunwtd.dfr$resid)
                thelayout <- c(3, 1)
              }
              
              thetitle <- "Quantile-Quantile (QQ) Graphs on\nGaussian (Normal) Distribution"

              all.dfr$typef <- factorInSeq(all.dfr$type)

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
                                   top=list(pad1=0.5, pad2=1),
                                   right=list(pad1=0.5, pad2=2)
                                   ))
              on.exit(trellis.par.set("axis.components", trellisparstg3),
                      add=TRUE)              

              thegraph <- qqmath(~ resid | typef,
                                 distribution=qnorm,
                                 data=all.dfr,
                                 prepanel=prepanel.qqmathline,
                                 panel=function(x, ...) {
                                   panel.grid(h=-1, v=-1)
				   panel.number <- panel.number()
                                   if(is.element(panel.number,
                                                 c(1, nlevels(all.dfr$typef)))) {
                                     panelside <- ifelse(panel.number==1,
                                                         "left", "right") 
                                     tickmarks <- setupAxisTicks(all.dfr$resid,
                                                                 logscale=FALSE,
                                                                 grid=TRUE)
                                     panel.axis(side=panelside,
                                                at=tickmarks,
                                                labels=names(tickmarks),
                                                tck=0.20, text.cex=0.6,
                                                rot=0,
                                                outside=TRUE)
                                   }
                                   snq <- qnorm(ppoints(length(x)))
                                   themin <- max(snq) - 1.04*(max(snq) - min(snq))
                                   panel.text(x=rep(themin, 2),
                                              y=range(x),
                                              labels=paste("",
                                                c(signif(min(x), 3),
                                                  signif(max(x), 3))),
                                              col="blue", adj=0, cex=0.6)

                                   if(is.null(line) || line) {
                                     panel.qqmathline(x, ...)
                                   }
                                   panel.qqmath(x, ...)
                                 },
                                 col="black",
                                 xlab="Standard Gaussian (Normal) Quantile",
                                 ylab=list(cex=0.7,
                                   label="  "),
                                 layout=thelayout, aspect=1,
                                 scales=list(y=list(labels=NULL, tck=0),
                                   x=list(tck=0.20, cex=0.7)),
                                 main=list(label=paste(thetitle, "\n",
                                             analysisname, sep=""), cex=1.1),
                                 par.strip.text=list(cex=0.7))
              print(thegraph) ## , position=c(0,0,0.95,1))

              upViewport(0)
              seekViewport(trellis.vpname("ylab"))
              ylabchar <- catCharExpr("Residual in Fit of ", endptlabel)
              grid.text(ylabchar,
                        x = unit(0, "lines"), rot=90, gp=gpar(cex=0.9))
              upViewport(0)
              
              if(settings$stamps) graphStampCG()
            }

            else if(model=="olsonly" && ols && device=="single") {
              qqgraph(resid(olsfit),
                      analysisname = analysisname,
                      desc="Classical fit",
                      endptname=catCharExpr(" Fit of ", endptlabel),
                      line=line,
                      titlestamp = TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)
            }

            else if((model=="rronly" || model=="rrwtdonly") && rr
                    && device=="single") {
              qqgraph(sqrt(rrfit$w)*resid(rrfit),
                      line=line,
                      desc = "Resistant & Robust fit Weighted",
                      analysisname = analysisname, 
                      endptname=catCharExpr(" Fit of ", endptlabel),
                      titlestamp = TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)
            }

            else if(model=="rrunwtdonly" && rr && device=="single") {
              qqgraph(resid(rrfit),
                      line=line,
                      desc = "Resistant & Robust fit Unweighted",
                      analysisname = analysisname, 
                      endptname=catCharExpr(" Fit of ", endptlabel),
                      titlestamp = TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)
            }

            else if(rr && ols &&
                    is.element(device, c("ask","multiple")) &&
                    is.element(model, c("both", "extended"))) {
              if(device=="ask") {
                op <- par(ask = TRUE)
                on.exit(par(op), add=TRUE)
              }
              else { ## device=="multiple"
                cat(cgMessage("New graphic devices will be generated to",
                              "hold multiple graphs requested...\n",
                              warning=TRUE))
              }
              qqgraph(resid(olsfit),
                      line=line,
                      desc = "Classical",
                      analysisname = analysisname, 
                      endptname=catCharExpr(" Fit of ", endptlabel),
                      titlestamp = TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)
              
              if(device=="multiple") {
                do.call("cgDevice", c(list(new=TRUE), dots))
              }
              qqgraph(sqrt(rrfit$w)*resid(rrfit),
                      line=line,
                      desc = "Resistant & Robust fit Weighted",
                      analysisname = analysisname, 
                      endptname=catCharExpr(" Fit of ", endptlabel),
                      titlestamp = TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)

              if(model=="extended") {
                if(device=="multiple") {
                  do.call("cgDevice", c(list(new=TRUE), dots))
                }
                qqgraph(resid(rrfit),
                        line=line,
                        desc = "Resistant & Robust fit Unweighted",
                        analysisname = analysisname, 
                        endptname=catCharExpr(" Fit of ", endptlabel),
                        titlestamp = TRUE)
                if(settings$stamps) graphStampCG(grid=FALSE)
              }
            }            

            else if(uv && device=="single") {
              qqgraph(residuals(uvfit, type="pearson"),
                      line=line,
                      desc = "Unequal Variances fit",
                      analysisname = analysisname, 
                      endptname=catCharExpr("  , Standardized  ",
                        catCharExpr(" Fit of ", endptlabel), rev=TRUE),
                      titlestamp = TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)
            }

            else if(aft && device=="single") {
              qqgraph(residuals(aftfit, type="response"),
                      status=aftfit$dfru$status,
                      line=line,
                      desc = "AFT Censored",
                      analysisname = analysisname, 
                      endptname=catCharExpr(" Fit of ", endptlabel),
                      titlestamp = TRUE)
              if(settings$stamps) graphStampCG(grid=FALSE)
            }

            else {
              stop(cgMessage("The chosen device and model arguments",
                             "are not compatible either with each other",
                             "or with the fitted model(s) in the fit object.",
                             seeHelpFile("qqGraph")))
            }

            invisible()
          })

setClass("cgOneFactorDownweightedTable",
         representation(contents="dataframeOrNULL", 
                        settings="list"),
         prototype(contents=NULL, settings=list()))

setMethod("downweightedTable", "cgOneFactorFit",
          downweightedTable.cgOneFactorFit <- 
          function(fit, cutoffwt, display="print", ...) {
            ##
            ## PURPOSE: Create a table of data observations
            ## that were downweighted by MASS::rlm()
            ##
            ## Input arguments check
            if(class(rrfit <- fit@rrfit)[1]!="rlm") {
              stop(cgMessage("No resistant & robust fit of class rlm",
                             "is available."))
            }
            
            validCutoffWt(cutoffwt)
            display <- validArgMatch(display, c("print","none","show"))

            dots <- list(...)
            validDotsArgs(dots, names=c(""))
            ## End input arguments check

            settings <- fit@settings
            settings$cutoffwt <- cutoffwt
            tbl <- message <- NULL ## intializations
            
            obsweights <- sqrt(rrfit$w)
            flaggedpoints <- (obsweights <= cutoffwt)

            if(any(flaggedpoints)) {
              tbl <- data.frame(rrfit$dfru[flaggedpoints, ],
                                obsweights[flaggedpoints],
                                100*(1 - obsweights[flaggedpoints]),
                                row.names=NULL)
              
              names(tbl) <- c("group", "endpoint", "weight", "pct down-weighted")
            }

            returnObj <- new("cgOneFactorDownweightedTable",
                             contents=tbl, 
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


setMethod("print", "cgOneFactorDownweightedTable",
          print.cgOneFactorDownweightedTable <-
          function(x, digits=NULL, title=NULL, endptname=NULL, ...) {
            ##
            ## PURPOSE: Semi-formatted print version of Downweighted Table
            ## from Resistant & Robust fit
            ##
            ## NOTE: Had to use x as an argument because of the system defined
            ## generic. I would have preferred to use object; hence the first
            ## statement below.
            ##
            object <- x
            settings <- object@settings
            table <- object@contents
            cutoffwt <- settings$cutoffwt

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
              title <- paste("Downweighted Observations Table from ",
                             "Resistant & Robust Fit\n",
                             settings$analysisname, sep="") 
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

            if(!is.null(table)) {
              table$group <-   format(table$group, justify="left")
              table$endpoint <- fround(table$endpoint, digits)
              table$weight <- fround(table$weight, digits=2)
              table$"pct down-weighted" <- fmtPercent(table$"pct down-weighted")
              
              curwidth <- getOption("width")
              on.exit(options(width=curwidth), add=TRUE)
              if(curwidth < 500) { options(width=500) }

              cat(paste("\nWeights less than ", cutoffwt,
                        "\n(i.e. Observations Downweighted by ", fmtPercent(100*(1-cutoffwt)),
                        "% or more", ")",
                        "\n\n",
                        sep=""))
              print(table, quote=FALSE, row.names=FALSE)
            }
            else cat("\n","No observations have weights less than\n",
                     "the cutoff of ",
                     cutoffwt, ", so the table is empty.",
                     "\n\n(i.e. ", fmtPercent(100*(1-cutoffwt)),
                     "% or more of a downweight of an observation did not occur.", ")",
                     "\n",
                     sep="")
            
            invisible()
          })

setMethod("show", "cgOneFactorDownweightedTable",
          show.cgOneFactorDownweightedTable <- function(object) showDefault(object))


validCutoffWt <- function(x) {
  if(missing(x) ||  !(x > 0 && x < 1) ) {
    stop(cgMessage("The cutoffwt value for identifying",
                   "downweighted points",
                   "needs to be numerically specified",
                   "and also be greater than zero and less than",
                   "one.",
                   seeHelpFile("downweightedTable")))
  }
  
  return(TRUE)
}









