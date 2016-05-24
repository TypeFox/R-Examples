## $Id:  $
## Paired Groups Case

## Fit Paired Groups Data

setClass("cgPairedDifferenceFit",
         representation(olsfit="olsfit",
                        rrfit="rrfit",                 
                        settings="list"))

setMethod("fit", "cgPairedDifferenceData",
          fit.cgPairedDifferenceData <- function(data, type="rr", ...) {
            ##
            ## PURPOSE: Fit data that has a 2-group paired structure
            ##
            ## NOTES: Paired Differences are modeled directly.
            ##
            ## Input arguments handling
            dots <- list(...)
            validDotsArgs(dots, names=c("maxIter"))
            
            ## Initializations ('ols' is always TRUE so not specified)
            rr <- FALSE
            rrfit <- "No fit was requested."

            if(missing(type)) type <- "rr"  ## due to fit generic definition
            
            if(class(try(match.arg(type, c("ols","rr"))))=="try-error") {
              stop(cgMessage("The type argument needs to evaluate to one",
                             "of: \"ols\", \"rr\".",
                             seeHelpFile("CGPairedDifferenceFit")))
            }
            
            maxIter <- if(is.null(dots$maxIter)) 100 else dots$maxIter 
            validNumeric(maxIter, positive=TRUE, integer=TRUE)

            if(type=="rr") { rr <- TRUE }

            ## Need to check about type="ols"
            ## End input arguments handling

            dfr.gcfmt <- data@dfr.gcfmt
            dfru<- data@dfru
            
            settings <- data@settings  
            endptscale <- settings$endptscale
            grpnames <- settings$grpnames

            oldop <- options(contrasts=c("contr.treatment", "contr.poly"))
            on.exit(oldop, add=TRUE)

            endpt <- with(dfr.gcfmt,
                          if(endptscale=="log") {
                            difflogendpt
                          }
                          else {
                            diffendpt
                          })
                          
            ## Ordinary Least Squares is *always* fit
            olsfit <- lm(endpt ~ 1)
            olsfit$dfr.gcfmt <- dfr.gcfmt
            olsfit$dfru <- dfru

            ## May also need a two-sample fit that ignores the pairing for use by
            ## samplesize methods later
            ## Note that "endpt" refers to "dfru$endpt" in these lm calls.
            olsfit$onefactorolsfit <- if(endptscale=="log") {
              lm(log(endpt) ~ -1 + grpf, data=dfru)
            }
            else {
              lm(endpt ~ -1 + grpf, data=dfru)
            }

            ## lme version to help assess between and within variance
            ## components.
            olsfit$lmefit <- if(endptscale=="log") {
              lme(log(endpt) ~ grpf, data=dfru,
                  random=list(expunitf= ~ 1), method="REML")
            }
            else {
              lme(endpt ~ grpf, data=dfru,
                  random=list(expunitf= ~ 1), method="REML")
            }

            ## Thanks to tip from Gang Chen at NIH on R-help
            ## https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q2/016300.html
            ## Gets individual estimates
            olsfit$lmefitindv <- if(endptscale=="log") {
              lme(log(endpt) ~ -1 + grpf, data=dfru,
                  random=list(expunitf=pdCompSymm(~ -1 + grpf)), method="REML")
            }
            else {
              lme(endpt ~ -1 + grpf, data=dfru,
                 random=list(expunitf=pdCompSymm(~ -1 + grpf)), method="REML")
            }
                        
            if(rr) {  
              rrfit <- try(rlm(endpt  ~ 1, method="MM", maxit=maxIter))
              if(class(rrfit)[1]!="try-error") {
                ## that is, convergence occurred
                rrfit$dfr.gcfmt <- dfr.gcfmt
                rrfit$dfru <- dfru
              
                rrfit$onefactorfit <- if(endptscale=="log") {
                  try(rlm(log(endpt) ~ -1 + grpf, data=dfru,
                          method="MM", maxit=maxIter))
                }
                else {
                  try(rlm(endpt ~ -1 + grpf, data=dfru,
                          method="MM", maxit=maxIter))
                }

                ## lme version to help assess between and within variance
                ## components, using weighted versions of data based
                ## on rrfit$w weights
                rrfit$lmefit <- if(endptscale=="log") {
                  lme(log(endpt) ~ -1 + grpf,
                      data=data.frame(dfru, w=rep(rrfit$w, 2)),
                      random=list(expunitf=pdSymm(~ -1 + w)), method="REML")
                }
                else {
                  lme(endpt ~ -1 + grpf,
                      data=data.frame(dfru, w=rep(rrfit$w, 2)),
                      random=list(expunitf=pdSymm(~ -1 + w)), method="REML")
                }
                
                ## workaround so return object instance below will
                ## accept rrfit in the slot. Has something to do with
                ## setClassUnion("rrfit", c("character", "rlm", "lm")) idiom I use.
                ## see https://stat.ethz.ch/pipermail/r-devel/2008-April/049278.html
                class(rrfit) <- "rlm"
                
              }
              else {
                warning(cgMessage("The Resistant & Robust (rr) fit did not",
                                  "converge in the specified number of",
                                  "iterations. You may want to try again with an",
                                  "increased value for the maxIter argument.",
                                  warning=TRUE))
                rrfit <- "The fit did not converge in the specified number of iterations."
              }
            }
            
            returnObj <- new("cgPairedDifferenceFit",
                             olsfit=olsfit,
                             rrfit=rrfit,
                             settings=settings)

            returnObj
          })


setMethod("print", "cgPairedDifferenceFit",
          print.cgPairedDifferenceFit <-
          function(x, title=NULL, endptname=NULL,...) {
            ##
            ## PURPOSE: Simple print version of cgPairedDifferenceFit
            ## object. Echoes print methods for individual
            ## object classes, such as lm and rlm
            ##
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
            
            olsfit <- x@olsfit
            rrfit <- x@rrfit

            settings <- x@settings

            ols <- rr <- FALSE  ## Initializations
            if(class(rrfit)[1]=="rlm" && model!="olsonly") {
              rr <- TRUE
            }
            if(class(olsfit)[1]=="lm" && (model!="rronly")) {
              ols <- TRUE
              if(!rr) model <- "olsonly"
            }

            if(is.null(title)) {
              title <- paste("Fitted Models of", settings$analysisname) 
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
            if(endptname!="") {
              grpnames <- settings$grpnames
              refgrp <- settings$refgrp
              minuendgrp <- grpnames[grpnames != refgrp]
              cat(paste("Endpoint:",
                        minuendgrp, "vs.", refgrp,
                        "Difference in",
                        if(settings$endptscale=="log") "log",
                        endptname, "\n")) }
            if(ols) {
              cat("\nClassical Least Squares Model Fit\n")
              print(olsfit, ...)
            }
            
            if(rr) {
              cat("\nResistant & Robust Model Fit\n\n")
              print(rrfit, ...)
              cat("\nEstimated Standard Deviation from rlm is",
                  signif(summary(rrfit)$stddev, 4), "\n\n")              
            }

            invisible()  
          }
          )

setMethod("showObj", "cgPairedDifferenceFit",
          showObj.cgPairedDifferenceFit <- function(object) showDefault(object))

setMethod("summary", "cgPairedDifferenceFit",
          summary.cgPairedDifferenceFit <-
          function(object, title=NULL, endptname=NULL, ...) {
            ##
            ## PURPOSE: Simple summary of cgPairedDifference fit
            ## objects that are focused on the paired difference measures.
            ## Echoes summary methods for individual
            ## object classes lm and rlm
            ##
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
            
            olsfit <- object@olsfit
            rrfit <- object@rrfit

            settings <- object@settings

            ols <- rr <- FALSE  ## Initializations
            if(class(rrfit)[1]=="rlm" && model!="olsonly") {
              rr <- TRUE
            }
            if(class(olsfit)[1]=="lm" && (model!="rronly")) {
              ols <- TRUE
              if(!rr) model <- "olsonly"
            }

            if(is.null(title)) {
              title <- paste("Fitted Model Summaries of", settings$analysisname) 
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
            if(endptname!="") {
              grpnames <- settings$grpnames
              refgrp <- settings$refgrp
              minuendgrp <- grpnames[grpnames != refgrp]
              cat(paste("Endpoint:",
                        minuendgrp, "vs.", refgrp,
                        "Difference in",
                        if(settings$endptscale=="log") "log",
                        endptname, "\n")) }
            
            if(ols) {
              cat("\nClassical Least Squares Model Fit Summary\n")
              print(summary(olsfit, ...))
            }
            
            if(rr) {
              cat("\nResistant & Robust Model Fit Summary\n")
              print(summary(rrfit, ...))
              cat("\nEstimated Standard Deviation from rlm is",
                  signif(summary(rrfit)$stddev, 4), "\n\n")
            }

            invisible()  

          }
          )


