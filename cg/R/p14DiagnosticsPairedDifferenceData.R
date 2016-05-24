## $Id:  $
## Paired Groups

## Diagnostic methods

## There is no need for a varianceGraph method  because the formal analysis is
## essentially on one group of data, i.e. the paired differences.  

setMethod("qqGraph", "cgPairedDifferenceFit",
          qqGraph.cgPairedDifferenceFit <-
          function (fit, line = TRUE, cgtheme = TRUE, device = "single", ...)
          {
            ## PURPOSE: Create a Q-Q Gaussian Plot(s) of the Residuals
            ##
            ## NOTES: The definition of qqgraph() is in
            ## p04DiagnosticsOneFactorData.R
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names=c("model"))
            validBoolean(cgtheme)

            rrfit <- fit@rrfit
            olsfit <- fit@olsfit

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
            
            ols <- rr <- FALSE  ## Initializations
                        
            if(class(rrfit)[1]=="rlm" && model!="olsonly") {
              rr <- TRUE
            }
            if(class(olsfit)[1]=="lm" && model!="rronly") {
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
              print(thegraph) 

              upViewport(0)
              seekViewport(trellis.vpname("ylab"))
              ylabchar <- catCharExpr("Residual in Fit of ",
                                      endptlabel)
              
              grid.text(ylabchar,
                        x = unit(0, "lines"), rot=90, gp=gpar(cex=0.9))
              grid.text("(Paired Difference)",
                        x = unit(1, "lines"), rot=90, gp=gpar(cex=0.7))
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
              mtext(side=2, line=2.5, text="(Paired Difference)", cex=0.7)
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
              mtext(side=2, line=2.5, text="(Paired Difference)", cex=0.7)
              if(settings$stamps) graphStampCG(grid=FALSE)
            }

            else if(model=="rrunwtdonly" && rr && device=="single") {
              qqgraph(resid(rrfit),
                      line=line,
                      desc = "Resistant & Robust fit Unweighted",
                      analysisname = analysisname, 
                      endptname=catCharExpr(" Fit of ", endptlabel),
                      titlestamp = TRUE)
              mtext(side=2, line=2.5, text="(Paired Difference)", cex=0.7)
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
              mtext(side=2, line=2.5, text="(Paired Difference)", cex=0.7)
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
              mtext(side=2, line=2.5, text="(Paired Difference)", cex=0.7)
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
                mtext(side=2, line=2.5, text="(Paired Difference)", cex=0.7)
                if(settings$stamps) graphStampCG(grid=FALSE)
              }
            }            
            else {
              stop(cgMessage("The chosen device and model arguments",
                             "are not compatible either with each other",
                             "or with the fitted model(s) in the fit object.",
                             seeHelpFile("qqGraph")))
            }
                       
            invisible()
          }
          )

setClass("cgPairedDifferenceDownweightedTable",
         representation(contents="dataframeOrNULL",
                        settings="list"),
         prototype(contents=NULL, settings=list()))

setMethod("downweightedTable", "cgPairedDifferenceFit",
          downweightedTable.cgPairedDifferenceFit <-
          function (fit, cutoffwt, display = "print", ...) {
            ##
            ## PURPOSE:Create a table of paired difference observations
            ## that were downweighted by MASS::rlm()
            ##
            ## NOTES: The definition of validCutOffWt() is in
            ## p04DiagnosticsOneFactorData.R
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
            tbl <- message <- NULL
            
            obsweights <- sqrt(rrfit$w)
            flaggedpoints <- (obsweights <= cutoffwt)

            if(any(flaggedpoints)) {
              dfr.gcfmt <- rrfit$dfr.gcfmt              
              tbl <- data.frame(dfr.gcfmt[flaggedpoints, 1:3],
                                obsweights[flaggedpoints],
                                100*(1 - obsweights[flaggedpoints]),
                                dfr.gcfmt[flaggedpoints, 4],
                                row.names=NULL)

              if(settings$endptscale== "log") {
               tbl <- data.frame(tbl, cbind(exp(dfr.gcfmt[flaggedpoints, 5]),
                      100*(exp(dfr.gcfmt[flaggedpoints, 5]) - 1)))
              }

              names(tbl) <- c(settings$expunitname, settings$grpnames,
                              "weight", "pct down-weighted",
                              "simple diff",
                              if(settings$endptscale== "log") {
                                c("ratio diff", "pct diff") }
                              )
            }
            
            returnObj <- new("cgPairedDifferenceDownweightedTable",
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

setMethod("print", "cgPairedDifferenceDownweightedTable",
          print.cgPairedDifferenceDownweightedTable <-
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
              title <- paste("Downweighted Paired Differences Table from ",
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
              fmt.table <- table[, 1:3]
              fmt.table[,1] <- format(table[,1], justify="left")
              fmt.table[,2] <- fround(table[,2], digits)
              fmt.table[,3] <- fround(table[,3], digits)
              fmt.table$weight <- fround(table$weight, digits=2)
              fmt.table$"pct down-weighted" <- fmtPercent(table$"pct down-weighted")
              fmt.table$"Simple Diff" <- fround(table$"simple diff",
                                                      digits)
                           
              if(settings$endptscale == "log") {
                fmt.table$"Ratio Diff" <- fmtRatio(table$"ratio diff")
                fmt.table$"Pct Diff" <- fmtPercent(table$"pct diff")
              }

              curwidth <- getOption("width")
              on.exit(options(width=curwidth), add=TRUE)
              if(curwidth < 500) { options(width=500) }

              cat(paste("\nWeights less than ", cutoffwt,
                        "\n(i.e. Paired Differences Downweighted by ",
                        fmtPercent(100*(1-cutoffwt)),
                        "% or more", ")",
                        "\n\n",
                        sep=""))
              print(fmt.table, quote=FALSE, row.names=FALSE)
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

setMethod("show", "cgPairedDifferenceDownweightedTable",
          show.cgPairedDifferenceDownweightedTable <- function(object) showDefault(object))

