## $Id:  $
## Paired Groups Case

setClass("cgPairedDifferenceComparisonsTable",
         representation(ols.comprs="dataframeMatrixOrNULL",
                        rr.comprs="dataframeMatrixOrNULL",
                        settings="list"),
         prototype(ols.comprs=NULL,
                   rr.comprs=NULL,
                   settings=list()))

setMethod("comparisonsTable", "cgPairedDifferenceFit",
          comparisonsTable.cgPairedDifferenceFit <-
          function(fit, type="pairwisereflect",
                   alpha=0.05, addpct=FALSE, 
                   display="print", ...) {
            ##
            ## PURPOSE: Construct comparison table of info for single paired difference
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names="model")
            type <- validArgMatch(type, c("pairwisereflect","pairwise"))

            validAlpha(alpha)
            validBoolean(addpct)

            display <- validArgMatch(display, c("print","none","show"))

            modelarg <- getDotsArgName(dots, "model")
            if(!is.na(modelarg)) {
              model <- eval(parse(text=paste("dots$", modelarg, sep="")))
              model <- validArgMatch(model, choices=c("both", "olsonly","rronly"))
            }
            else {
              model <- "both"
            }

            ## The following arguments defined in the generic are not needed
            ## and thus ignored
            mcadjustarg <- getDotsArgName(dots, name="mcadjust")
            if(!is.na(mcadjustarg)) {
              cat(cgMessage("The mcadjust argument seems to be specified,",
                            "perhaps by an unnamed argument value. It is not applicable for this",
                            "method and is ignored...\n\n", warning=TRUE))
            }
            contrastmatrixarg <- getDotsArgName(dots, name="contrastmatrix")
            if(!is.na(contrastmatrixarg)) {
              cat(cgMessage("The contrastmatrix seems to be specified,",
                            "perhaps by an unnamed argument value. It is not applicable for this",
                            "method and is ignored...\n\n", warning=TRUE))
            }
            refgrparg <- getDotsArgName(dots, name="refgrp")
            if(!is.na(refgrparg)) {
              cat(cgMessage("The refgrp argument seems to be specified,",
                            "perhaps by an unnamed argument value. It is not applicable for this",
                            "method and is ignored...\n\n", warning=TRUE))
            }            

            ## initializations
            ols <- rr <- FALSE
            ols.comprs <- rr.comprs <- NULL
            ##
            settings <- fit@settings
            endptscale <- settings$endptscale
            rrfit <- fit@rrfit
            olsfit <- fit@olsfit
            grpnames <- settings$grpnames
            offset <- settings$addconstant

            validAddPct(addpct, endptscale)

            if(class(rrfit)[1]=="rlm" && model!="olsonly") {
              rr <- TRUE
            }
            if(class(olsfit)[1]=="lm" && model!="rronly") {
              ols <- TRUE
              if(!rr) model <- "olsonly"
            }

            df.residual <- olsfit$df.residual
            
            ## proof of concept function for comparisons with essentially
            ## just one group, the paired differences. But we also need to
            ## take into account that we can get indvidual group means
            ## for the Classical least squares case through the alternative
            ## lme fit. No such analogue exists for the resistant/robust case
            ## Adapt the comparisons function, and long-term try to merge this into
            ## it
            comparisons.paired <- function(diffestimate,
                                           var.diffest,
                                           model="ols",
                                           estimates,
                                           varcovmatrix, errordf=Inf,
                                           endptscale, alpha=0.05, 
                                           type="pairwisereflect",
                                           cnames,
                                           analysisname="", endptname="",
                                           digits=NULL, addpct=FALSE,
                                           display="print", ...) {
              ## Input argument checking
              validAlpha(alpha)
              validBoolean(addpct)
              type <- validArgMatch(type, c("pairwisereflect","pairwise"))
              model <- validArgMatch(model, c("ols","rr"))
              
              if(model=="ols") {
                if(!hasArg(estimates)) reportInvalidArg("estimates")
                if(!hasArg(varcovmatrix)) reportInvalidArg("varcovmatrix")              
                if(missing(cnames) || !(length(cnames)==2 && validCharacter(cnames))) {
                  stop(cgMessage("The cnames argument needs to be of",
                                 "length 2 and be of character type."))
                  estimates <- validEstimates(estimates)
                }
              }

              validNumeric(diffestimate)
              validNumeric(var.diffest)
              
              endptscale <- validArgMatch(endptscale, c("log","original"))
              validAddPct(addpct, endptscale)
              
              display <- validArgMatch(display, c("print","none","show"))

              if(type=="pairwisereflect") {
                dL <- matrix(c(1, -1), nrow=2)
                rownames(dL) <- c(paste(cnames[2], "vs.", cnames[1]),
                                  paste(cnames[1], "vs.", cnames[2]))
              }
              else { ## type=="pairwise"
                dL <- matrix(1) ## 1x1 matrix
                rownames(dL) <- paste(cnames[2], "vs.", cnames[1])
              }

              estdiff <- dL %*% matrix(diffestimate)
              sediff <- sqrt(diag(dL %*%  matrix(var.diffest)  %*% t(dL)))

              tcrit <- qt(1 - alpha/2, errordf)
              lowerci <- estdiff - tcrit * sediff
              upperci <- estdiff + tcrit * sediff
              pval <- 2 * (1 - pt(abs(estdiff/sediff), errordf))

              mcp <- data.frame(estimate=estdiff, se=sediff,
                                lowerci=lowerci, upperci=upperci,
                                pval=pval)

              if(endptscale=="log") {
                logscalest <- mcp[,1]
                mcp[, 1] <- 100 * ( exp(logscalest) - 1 ) # Point Estimate
                mcp[, 2] <- 100 * ( exp(logscalest) * mcp[, 2] ) # Std Err
                mcp[, 3:4] <- 100 * ( exp(mcp[, 3:4]) - 1) # Confidence Limits
              }

              ## Add Individual group means and standard errors if ols
              if(model=="ols") {
                seindv <- sqrt(diag(varcovmatrix))

                if(nrow(mcp)==2) {
                  mcp$meanA <- rev(estimates)
                  mcp$seA <- seindv
                  mcp$meanB <- estimates
                  mcp$seB <- seindv
                }
                else { ## nrow(mcp)==1
                  mcp$meanA <- estimates[2]
                  mcp$seA <- seindv[2]
                  mcp$meanB <- estimates[1]
                  mcp$seB <- seindv[1]                  
                }

                if(endptscale=="log") {
                  ## Now focus on the individual estimates
                  if(is.null(offset)) {
                    mcp$meanA <- exp(mcp$meanA)
                    mcp$meanB <- exp(mcp$meanB)
                  }
                  else { ## offset
                    mcp$meanA <- exp(mcp$meanA) - offset
                    mcp$meanB <- exp(mcp$meanB) - offset
                  }
                  
                  mcp$seA <- mcp$meanA * mcp$seA
                  mcp$seB <- mcp$meanB * mcp$seB
                }

                ## and if original scale is used, calculate
                ## percent difference(s) from the arithmetic means
                if(endptscale=="original" && addpct==TRUE) {
                  mcp$pctdiff <- 100*(mcp$meanA / mcp$meanB - 1)
                }

              }

              if(display=="print") {
                fmt.mcp <- mcp
                cat(paste("Comparisons Table",
                          if(analysisname!="") paste("for", analysisname), "\n"))
                if(endptname!="") { cat(paste("Endpoint:", endptname, "\n")) }
                diffmetric <- "Differences (A vs. B)"
                if(endptscale=="log") {
                  diffmetric <- paste("Percent", diffmetric)
                }
                cat(diffmetric, "\n")

                pvalfmt <- fmtPvalue(fmt.mcp$pval)
                if(regexpr("Percent", diffmetric) > 0) {
                  fmt.mcp$estimate <- fmtPercent(fmt.mcp$estimate)
                  fmt.mcp$se <- fmtPercent(fmt.mcp$se)
                  fmt.mcp$lowerci <- fmtPercent(fmt.mcp$lowerci)
                  fmt.mcp$upperci <- fmtPercent(fmt.mcp$upperci)
                  fmt.mcp$pval <- pvalfmt
                  if(!is.null(digits) && model=="ols" ) {
                    fmt.mcp$meanA <- fround(fmt.mcp$meanA, digits=digits)
                    fmt.mcp$seA <- fround(fmt.mcp$seA, digits=digits)
                    fmt.mcp$meanB <- fround(fmt.mcp$meanB, digits=digits)
                    fmt.mcp$seB <- fround(fmt.mcp$seB, digits=digits)
                    
                    { names(fmt.mcp)[is.element(names(fmt.mcp),
                                                c("meanA","seA",
                                                  "meanB","seB"))] <- c("geomeanA","seA",
                                                                        "geomeanB","seB") }
                    { names(mcp)[is.element(names(mcp), c("meanA","seA",
                                                          "meanB","seB"))] <-
                                                            c("geomeanA","seA",
                                                              "geomeanB","seB") }
                  }
                }
                else { ## Simple Differences
                  pctdiff <- mcp$pctdiff
                  if(!is.null(digits)) fmt.mcp <- fround(fmt.mcp, digits)
                  fmt.mcp$pval <- pvalfmt
                  if(endptscale=="original" && addpct==TRUE && model=="ols") {
                    fmt.mcp$pctdiff <- fmtPercent(pctdiff)
                  }
                }

                cat(paste(round(100*(1-alpha), 0), "% Confidence ",
                          "(alpha of ", alpha, ")",
                          "\n",
                          sep=""))

                curwidth <- getOption("width")
                on.exit(options(width=curwidth), add=TRUE)
                if(curwidth < 500) { options(width=500) }
                
                print(fmt.mcp, quote=FALSE)
              }
              else if (display=="show"){
                showDefault(mcp)
              }
              ## else show nothing
              ## return
              invisible(mcp)
            }

            if(ols) {
              olsdiffestimate <- olsfit$coef
              varcov.olsdiffestimate <- vcov(olsfit)
              olsindvestimates <- fixef(olsfit$lmefitindv)
              varcov.olsindvestimates <- vcov(olsfit$lmefitindv)

              ols.comprs <- comparisons.paired(olsdiffestimate,
                                               varcov.olsdiffestimate,
                                               model="ols",
                                               olsindvestimates,
                                               varcov.olsindvestimates,
                                               errordf=df.residual,
                                               endptscale=endptscale,
                                               alpha=alpha,
                                               type=type,
                                               cnames=grpnames,
                                               offset=offset,
                                               addpct=addpct, display="none")
            }
            if(rr) {
              rrdiffestimate <- rrfit$coef
              summ.rrfit <- summary(rrfit, method="XtWX", ...)
              stddev <- summ.rrfit$stddev
              varcov.rrdiffestimate <- (stddev^2) * summ.rrfit$cov.unscaled

              rr.comprs <- comparisons.paired(rrdiffestimate,
                                              varcov.rrdiffestimate,
                                              model="rr",
                                              errordf=df.residual,
                                              endptscale=endptscale,
                                              alpha=alpha,
                                              type=type,
                                              cnames=grpnames,
                                              addpct=addpct, display="none")
            }
            
            ## slightly different from the original fit object settings slot
            settings <- list(endptscale=settings$endptscale,
                             analysisname=settings$analysisname,
                             endptname=settings$endptname,
                             endptlabel=makeEndptLabel(settings$endptname,
                               settings$endptunits), 
                             alpha=alpha,
                             digits=settings$digits,
                             stamps=settings$stamps,
                             type=switch(type,
                               pairwise="Pairwise Comparisons (halfset)",
                               pairwisereflect="All Pairwise Comparisons"),
                             addpct=addpct)

            returnObj <- new("cgPairedDifferenceComparisonsTable",
                             ols.comprs=ols.comprs,
                             rr.comprs=rr.comprs,
                             settings=settings)
            
            if(display=="print") {
              print(returnObj)
            }
            else if(display=="show") {
              show(returnObj)
            }
            ## else display=="none"
            ## return
            invisible(returnObj)
            
          }
          )

setMethod("print", "cgPairedDifferenceComparisonsTable",
          print.cgPairedDifferenceComparisonsTable <-
          function (x, digits = NULL, title = NULL, endptname = NULL, ...) {
            ##
            ## PURPOSE: Semi-formatted print version of Comparisons Table
            ## 
            ## NOTE: Had to use x as an argument because of the system defined
            ## generic. I would have preferred to use object; hence the first
            ## statement below.
            ##
            object <- x
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

            olscomprs <- object@ols.comprs
            rrcomprs <- object@rr.comprs

            ols <- rr <- FALSE  ## Initializations
            if(!is.null(rrcomprs) && model!="olsonly") {
              rr <- TRUE
            }
            if(!is.null(olscomprs) && (model!="rronly")) {
              ols <- TRUE
            }

            settings <- object@settings
            alpha <- settings$alpha
            addpct <- settings$addpct
            
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
              title <- paste("Comparisons Table of", settings$analysisname) 
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

            diffmetric <- "Differences (A vs. B)"
            if(settings$endptscale=="log") {
              diffmetric <- paste("Percent", diffmetric)
            }
            cat(diffmetric, "\n")

            fmtdig <- function(x, diffmetric, digits, indvmeans=TRUE, addpct) {
              pvalfmt <- fmtPvalue(x$pval)
              fmt.x <- x
              if(regexpr("Percent", diffmetric) > 0) {
                fmt.x$estimate <- fmtPercent(x$estimate)
                fmt.x$se <- fmtPercent(x$se)
                fmt.x$lowerci <- fmtPercent(x$lowerci)
                fmt.x$upperci <- fmtPercent(x$upperci)
                fmt.x$pval <- pvalfmt
                if(indvmeans) {
                  fmt.x$meanA <- fround(x$meanA, digits=digits)
                  fmt.x$seA <- fround(x$seA, digits=digits)
                  fmt.x$meanB <- fround(x$meanB, digits=digits)
                  fmt.x$seB <- fround(x$seB, digits=digits)
                  { names(fmt.x)[is.element(names(fmt.x), c("meanA","seA",
                                                            "meanB","seB"))] <-
                                                              c("geomeanA","seA",
                                                                "geomeanB","seB") }
                }
              }
              else { ## Simple Differences
                fmt.x <- fround(x, digits)
                if(addpct==TRUE) {
                  pctdiff.fmt <- fmtPercent(x$pctdiff)
                  fmt.x$pctdiff <- pctdiff.fmt
                }
                fmt.x$pval <- pvalfmt
              }
              return(fmt.x)
            }

            informConfidence <- function() {
              cat(paste(round(100*(1-alpha), 0), "% Confidence ",
                        "(alpha of ", alpha, ")",
                        "\n",
                        sep=""))
            }

            curwidth <- getOption("width")
            on.exit(options(width=curwidth), add=TRUE)
            if(curwidth < 500) { options(width=500) }

            if(ols) {
              cat("\nClassical Least Squares Model Fit\n")
              informConfidence()
              print(fmtdig(olscomprs, diffmetric, digits, addpct=addpct), quote=FALSE)
            }
            
            if(rr) {
              cat("\nResistant & Robust Model Fit\n")
              informConfidence()
              print(fmtdig(rrcomprs, diffmetric, digits, indvmeans=FALSE,
                           addpct=FALSE)
                    , quote=FALSE)
            }
            
            invisible()
          })


setMethod("show", "cgPairedDifferenceComparisonsTable",
          show.cgPairedDifferenceComparisonsTable <- function(object)
          showDefault(object))


setMethod("comparisonsGraph", "cgPairedDifferenceComparisonsTable",
          comparisonsGraph.cgPairedDifferenceComparisonsTable <-
          function (compstable, cgtheme = TRUE, device = "single", wraplength = 20,
                    cex.comps = 0.7, ...) {
            ##
            ## PURPOSE: create a graph of CI's on paired difference
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names = c("model", "ticklabels"))
            validBoolean(cgtheme)

            settings <- compstable@settings
            difftype <- if (settings$endptscale == "log") "percent" else "simple"
            alpha <- settings$alpha
            alphapercent <- round(100 * alpha, 0)
            analysisname <- settings$analysisname
            endptlabel <- settings$endptlabel
            stamps <- settings$stamps
            digits <- settings$digits

            rrcomprs <- compstable@rr.comprs
            olscomprs <- compstable@ols.comprs[, 1:5]

            ticklabelsarg <- getDotsArgName(dots, "ticklabels")
            if(!is.na(ticklabelsarg)) {
              ticklabels <- eval(parse(text=paste("dots$", ticklabelsarg, sep="")))
              validList(ticklabels, names=c("mod","marks"),
                        argname="ticklabels")
            }
            else {
              ticklabels <- NULL
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

            ols <- rr <- FALSE  ## Initializations
            if(!is.null(rrcomprs) && model!="olsonly") {
              rr <- TRUE
            }
            if(!is.null(olscomprs) && (model!="rronly")) {
              ols <- TRUE
            }
            
            device <- validArgMatch(device, c("single","multiple", "ask"))

            thetitle <- "Comparisons Graph"
            if(rr && ols && is.element(model, "both") && device=="single") {
              all.dfr <- rbind(olscomprs, rrcomprs)
              all.dfr$typef <- factorInSeq(c(rep("Classical",
                                                 nrow(olscomprs)),
                                             rep("Resistant & Robust",
                                                 nrow(rrcomprs))))
              thetitle <- paste(thetitle, "s", sep="")
              numberofcomprs <- nrow(olscomprs)
              comprnames <- row.names(olscomprs)
              
              cgDevice(cgtheme=cgtheme)
              trellispanelstg <- trellis.par.get("clip")$panel
              trellis.par.set("clip", list(panel="off"))
              on.exit(trellis.par.set("clip", list(panel=trellispanelstg)),
                      add=TRUE)
              trellisparstg2 <- trellis.par.get("layout.widths")
              ymargin <- max(nchar(comprnames))/2
              trellis.par.set("layout.widths", list(left.padding=ymargin,
                                                    ylab=ymargin))
              on.exit(trellis.par.set("layout.widths", trellisparstg2),
                      add=TRUE)
              trellisparstg3 <- trellis.par.get("axis.components")
              trellis.par.set("axis.components",
                              list(left=list(pad1=0.2, pad2=2),
                                   bottom=list(pad1=0.2, pad2=1),
                                   top=list(pad1=1, pad2=1),
                                   right=list(pad1=0.5, pad2=1)
                                   ))
              on.exit(trellis.par.set("axis.components", trellisparstg3),
                      add=TRUE)


              allx <- unlist(all.dfr[, c("estimate","lowerci","upperci")])
              est <- all.dfr[, "estimate"]
              lower <- all.dfr[, "lowerci"]
              upper <- all.dfr[, "upperci"]
              compys <-  numberofcomprs:1

              if(difftype=="percent") {
                est <- log(pctToRatio(est))
                lower <- log(pctToRatio(lower))
                upper <- log(pctToRatio(upper))
                allx <- log(pctToRatio(allx))
                xmin <- min(c(0, allx))
                xmax <- max(c(0, allx))

                ## ideas in the next call adapted from errbar function
                thegraph <- xyplot(compys ~ Cbind(est,
                                                  lower,
                                                  upper) 
                                   | typef,
                                   data=all.dfr,
                                   digits=digits,
                                   panel = function(x, y, ...) {
                                     xlower <- attr(x, "other")[, 1]
                                     xupper <- attr(x, "other")[, 2]
                                     xcex <- 0.6
                                     xratioticks <-
                                       setupAxisTicks(exp(c(xlower, x, xupper)),
                                                      ratio=TRUE, percent=TRUE,
                                                      axis="x",
                                                      grid=TRUE, xcex=xcex,
                                                      digits=digits)
                                     if(!is.null(ticklabels)) {
                                       xratioticks <-
                                         makeTickMarks(ticklabels,
                                                       xratioticks,
                                                       percent=TRUE)
                                     }
                                     panel.axis(side="bottom",
                                                at=log(xratioticks),
                                                labels=names(xratioticks),
                                                text.cex=xcex,
                                                tck=0.20, outside=TRUE, rot=0)
                                     panel.points(x, compys, pch=16, col="black")
                                     lsegments(xlower, compys, xupper, compys)
                                     ##ycoord <- par()$usr[3:4]
                                     ##smidge <- 0.015 * (ycoord[2] - ycoord[1])
                                     smidge <- 0.015 * max(compys)
                                     lsegments(xlower, compys - smidge,
                                               xlower, compys + smidge)
                                     lsegments(xupper, compys - smidge,
                                               xupper, compys + smidge)
                                     ## zero reference line
                                     panel.abline(v=0, lty=3)
                                     
                                     ## Add 0 if not present
                                     if(all(xratioticks!=1)) {
                                       panel.axis(side="bottom",
                                                  at=0, labels="0", text.cex=0.6,
                                                  tck=0.20, outside=TRUE)
                                     }
                                     if(panel.number()==1) {
                                       ## Wrap long comparison labels
                                       for(i in seq(along=comprnames)) {
                                         if(nchar(comprnames[i]) > wraplength) {
                                           comprnames[i] <-
                                             paste(strwrap(comprnames[i],
                                                           width=wraplength,
                                                           exdent=1),
                                                   collapse="\n")
                                         }
                                       }
                                       panel.axis(side="left",
                                                  at=compys,
                                                  labels=comprnames,
                                                  tck=0.20, text.cex=cex.comps,
                                                  rot=0,
                                                  outside=TRUE)
                                     }
                                     xrange <- range(c(xlower, xupper))
                                     panel.text(y=rep(0.5, 2), x=xrange,
                                                labels=paste(
                                                  c(fmtRatioToPercent(exp(min(xrange))),
                                                    fmtRatioToPercent(exp(max(xrange)))), "\n"),
                                                col="blue", adj=0, cex=0.6)
                                   },
                                   layout=c(2,1), aspect=1, pch=16,
                                   col="black",
                                   as.table=TRUE,
                                   ylim = c(0.5, numberofcomprs + 0.5),
                                   xlim=rangeExtend(c(xmin, xmax), pct=10),
                                   xlab="Percent Difference",
                                   ylab="",
                                   scales=list(y=list(at=1:numberofcomprs,
                                                 labels=NULL, tck=c(0.15, 0), axs="i"),
                                     x=list(labels=NULL, tck=0)),
                                   main=list(label=paste(thetitle, "\n",
                                               analysisname, sep=""), cex=1.1),
                                   par.strip.text=list(cex=0.7)
                                   )
              }
              else {
                xmin <- min(c(0, allx))
                xmax <- max(c(0, allx))

                ## ideas in the next call adapted from Hmisc::errbar function
                thegraph <- xyplot(compys ~ Cbind(est,
                                                  lower,
                                                  upper) 
                                   | typef,
                                   data=all.dfr,
                                   digits=digits,
                                   panel = function(x, y, ...) {
                                     xlower <- attr(x, "other")[, 1]
                                     xupper <- attr(x, "other")[, 2]
                                     xcex <- 0.6
                                     xdiffticks <-
                                       setupAxisTicks(c(xlower, x, xupper),
                                                      difference=TRUE,
                                                      axis="x",
                                                      logscale=FALSE,
                                                      grid=TRUE,
                                                      xcex=0.6,
                                                      digits=digits)
                                     if(!is.null(ticklabels)) {
                                       xdiffticks <-
                                         makeTickMarks(ticklabels,
                                                       xdiffticks,
                                                       percent=FALSE)
                                     }
                                     panel.axis(side="bottom",
                                                at=xdiffticks,
                                                labels=names(xdiffticks),
                                                text.cex=xcex,
                                                tck=0.20, outside=TRUE, rot=0)
                                     panel.points(x, compys, pch=16, col="black")
                                     lsegments(xlower, compys, xupper, compys)
                                     ##ycoord <- par()$usr[3:4]
                                     ##smidge <- 0.015 * (ycoord[2] - ycoord[1])
                                     smidge <- 0.015 * max(compys)
                                     lsegments(xlower, compys - smidge,
                                               xlower, compys + smidge)
                                     lsegments(xupper, compys - smidge,
                                               xupper, compys + smidge)
                                     ## zero reference line
                                     panel.abline(v=0, lty=3)
                                     
                                     ## Add 0 if not present
                                     if(all(xdiffticks!=0)) {
                                       panel.axis(side="bottom",
                                                  at=0, labels="0", text.cex=0.6,
                                                  tck=0.20, outside=TRUE)
                                     }
                                     if(panel.number()==1) {
                                       ## Wrap long comparison labels
                                       for(i in seq(along=comprnames)) {
                                         if(nchar(comprnames[i]) > wraplength) {
                                           comprnames[i] <-
                                             paste(strwrap(comprnames[i],
                                                           width=wraplength,
                                                           exdent=1),
                                                   collapse="\n")
                                         }
                                       }
                                       panel.axis(side="left",
                                                  at=compys,
                                                  labels=comprnames,
                                                  tck=0.20, text.cex=cex.comps,
                                                  rot=0,
                                                  outside=TRUE)
                                     }
                                     xrange <- range(c(xlower, xupper))
                                     panel.text(y=rep(0.5, 2), x=xrange,
                                                labels=paste(
                                                  c(fmtDifference(min(xrange),
                                                                  digits=digits),
                                                    fmtDifference(max(xrange),
                                                                  digits=digits)), "\n"),
                                                col="blue", adj=0, cex=0.6)
                                   },
                                   layout=c(2,1), aspect=1, pch=16,
                                   col="black",
                                   as.table=TRUE,
                                   ylim=c(0.5, numberofcomprs + 0.5),
                                   xlim=rangeExtend(c(xmin, xmax), pct=10),
                                   xlab="Difference",
                                   ylab="",
                                   scales=list(y=list(at=1:numberofcomprs,
                                                 labels=NULL, tck=c(0.15, 0), axs="i"),
                                     x=list(labels=NULL, tck=0)),
                                   main=list(label=paste(thetitle, "\n",
                                               analysisname, sep=""), cex=1.1),
                                   par.strip.text=list(cex=0.7)
                                   )


              }
              
              print(thegraph) 
              seekViewport(trellis.vpname("xlab"))
              grid.text(catCharExpr("in ", endptlabel), y = unit(-1, "lines"))
              upViewport(0)

              comparisonsGraphStamp(mcadjust=FALSE, alphapercent, grid=TRUE,
                                    desc=settings$type)
              if(stamps) graphStampCG()

            }
            else if((model=="olsonly" || (ols && !rr && model=="both")) &&
                    device=="single") {
              comparisonsgraph(olscomprs, difftype, analysisname,
                               endptlabel, alpha, titlestamp=FALSE,
                               explanation=FALSE,
                               wraplength=wraplength,
                               cex.comps=cex.comps,
                               ticklabels=ticklabels)
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Comparisons Graph, Classical analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              comparisonsGraphStamp(mcadjust=FALSE, alphapercent,
                                    desc=settings$type)
            }
            
            else if((model=="rronly" || (!ols && rr && model=="both")) &&
                    device=="single") {
              ##            else if(model=="rronly" && rr && device=="single") {
              comparisonsgraph(rrcomprs, difftype, analysisname,
                               endptlabel, alpha, titlestamp=FALSE,
                               explanation=FALSE,
                               wraplength=wraplength,
                               cex.comps=cex.comps,
                               ticklabels=ticklabels)
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Comparisons Graph, Resistant & Robust analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              comparisonsGraphStamp(mcadjust=FALSE, alphapercent,
                                    desc=settings$type)
            }

            else if(rr && ols &&
                    is.element(device, c("ask","multiple")) &&
                    is.element(model, "both")) {
              
              device <- validArgMatch(device, c("multiple", "ask"))
              if(device=="ask") {
                op <- par(ask = TRUE)
                on.exit(par(op), add=TRUE)
              }
              
              comparisonsgraph(olscomprs, difftype, analysisname,
                               endptlabel, alpha, titlestamp=FALSE,
                               explanation=FALSE,
                               wraplength=wraplength,
                               cex.comps=cex.comps,
                               ticklabels=ticklabels)
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Comparisons Graph, Classical analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              comparisonsGraphStamp(mcadjust=FALSE, alphapercent,
                                    desc=settings$type)
              
              if(device=="multiple") {                
                eval(parse(text=paste("dots$", modelarg, " <- NULL", sep="")))
                do.call("cgDevice", c(list(new=TRUE), dots))
                cat(cgMessage("A new graphics device has been generated",
                              "to hold the Resistant & Robust",
                              "Comparisons Graph version.",
                              "The Classical Least Squares version is on the previous",
                              "device.\n",
                              warning=TRUE))
              }
              comparisonsgraph(rrcomprs, difftype, analysisname,
                               endptlabel, alpha, titlestamp=FALSE,
                               explanation=FALSE,
                               wraplength=wraplength,
                               cex.comps=cex.comps,
                               ticklabels=ticklabels)
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Comparison Graph, Resistant & Robust analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              comparisonsGraphStamp(mcadjust=FALSE, alphapercent,
                                    desc=settings$type)
            }

            else {
              stop(cgMessage("The chosen device and model arguments",
                             "are not compatible either with each other",
                             "or with the fitted model(s) in the fit object.",
                             seeHelpFile("comparisonsGraph")))
            }

            invisible()

          }
          )


setClass("cgPairedDifferenceVarianceTable",
         representation(contents="data.frame", efficiency="data.frame", settings="list"),
         prototype(contents=data.frame(), efficiency=data.frame(), settings=list()))


setMethod("varianceTable", "cgPairedDifferenceFit",
          varianceTable.cgPairedDifferenceFit <- function(fit,
                                                           display="print",
                                                           ...) {
            ##
            ## PURPOSE: Create a table of the estimated variance components
            ## (1) Within-subject variability due to pairings
            ## (2) Between-subject variability
            ## (Can only do the decomposition with Classical olsfit)
            ##
            display <- validArgMatch(display, c("print","none","show"))

            settings <- fit@settings
            expunitname <- settings$expunitname

            ## Extract the decomposition from the lmefit
            varcompest <-
              as.numeric(VarCorr(fit@olsfit$lmefit)[, 1])
            totalvar <- sum(varcompest)

            comps <- data.frame("Variance"=c(varcompest, totalvar),
                                "Percent"=100*c(varcompest/totalvar,1),
                                "Spread(StdDev)"=sqrt(c(varcompest,
                                  totalvar)),
                                check.names=FALSE)
            row.names(comps) <- c(paste("Within", expunitname),
                                  paste("Between", expunitname),
                                  "Total")

            ## Efficiency statement for using paired versus unpaired
            releff <- totalvar/varcompest[2]
            pctred.expunit <- 100*(1 - 1/releff)
            nexpunit <- nrow(fit@olsfit$dfr.gcfmt)
            releff.stmt <- paste("Assuming this data was collected under a paired",
                                 "design, the use of a Paired Difference analysis",
                                 "provides an estimated",
                                 paste(fmtPercent(pctred.expunit), "%", sep=""),
                                 "reduction in needed",
                                 paste(expunitname, "s", sep=""),
                                 "versus a unpaired",
                                 "(one factor) analysis. An estimated",
                                 floor(nexpunit*releff),
                                 paste(expunitname, "s", sep=""),
                                 "would have been needed under an unpaired design",
                                 "and analysis to achieve the same sensitivity as",
                                 "the",
                                 nexpunit,  paste(expunitname, "s", sep=""),
                                 "used in this Paired design and",
                                 "analysis.")

            efficiency <- data.frame(estimate=c(releff,pctred.expunit,
                                       nexpunit,floor(nexpunit*releff)))
            row.names(efficiency) <- c("Relative Efficiency",
                                       "Percent Reduction",
                                       paste("Number Paired ", expunitname, "s", sep=""),
                                       paste("Number Unpaired ", expunitname, "s",
                                             sep=""))
            attr(efficiency, "statement") <- releff.stmt
                                     
            returnObj <- new("cgPairedDifferenceVarianceTable",
                             contents=comps,
                             efficiency=efficiency,
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

setMethod("print", "cgPairedDifferenceVarianceTable",
          print.cgPairedDifferenceVarianceTable <-
          function(x, digits=NULL, title=NULL, endptname=NULL, ...) {
            ##
            ## PURPOSE: Semi-formatted print version of Variance Table
            ## 
            ## NOTE: Had to use x as an argument because of the system defined
            ## generic. Would have preferred to use object; hence the first
            ## statement below.
            ##
            object <- x
            settings <- object@settings

            if(is.null(digits)) {
              digits <- 2
            }
            curscipen <- getOption("scipen")
            on.exit(options(scipen=curscipen), add=TRUE)
            options(scipen=9)
            
            if(is.null(title)) {
              title <- paste("Variance Components Table of", settings$analysisname) 
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
            cat("\nClassical Least Squares Model Fit\n")

            out <- object@contents
            fmt.out <- out

            fmt.out[, c(1,3)] <- signif(out[, c(1,3)], digits=digits)
            fmt.out[, 2] <- fmtPercent(out[ , 2])
            
            curwidth <- getOption("width")
            on.exit(options(width=curwidth), add=TRUE)
            if(curwidth < 500) { options(width=500)}
            print(fmt.out, quote=FALSE)

            cat("\nNOTES:\n")
            
            cat(paragraphWrap(paste("1)", attr(object@efficiency, "statement"),"\n\n")))
            
            cat(paragraphWrap(paste("2) Analogous Resistant & Robust variance",
                                    "component estimates\n",
                                    "are not supported by the cg package.\n\n")))
                
            if(settings$endptscale=="log") {
              cat(paragraphWrap(paste("3) Since the fit is on log scale, the Spread(StdDev)",
                                      "column values can be\n",
                                      "interpreted as approximate",
                                      "Coefficients of Variation (CV) if\n",
                                      "less than 0.50.") ), "\n\n")
            }
            
            invisible()
          })

setMethod("show", "cgPairedDifferenceVarianceTable",
          show.cgPairedDifferenceVarianceTable <- function(object) showDefault(object))
