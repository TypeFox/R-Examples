## $Id: p01FitOneFactorData.R 3781 2013-01-11 20:07:34Z yye $
## One-Factor Unpaired Groups Case

## Fit One-Factor Unpaired Groups Data

setClass("cgOneFactorFit",
         representation(olsfit="olsfit",
                        rrfit="rrfit",
                        aftfit="aftfit",
                        uvfit="uvfit",
                        settings="list"))

setMethod("fit", "cgOneFactorData",
          fit.cgOneFactorData <- 
          function(data, type="rr", ...) {
            ##
            ## PURPOSE: Fit data that has a one-factor unpaired groups
            ## structure. Standard Least
            ## Squares Linear Model always fit, with a Resistant & Robust
            ## option based on M- and S-estimation. Also can fit
            ## an accelerated failure time model with lognormal or normal
            ## distribution on censored observations, or an unequal group
            ## variances structure.
            ## 
            ## Input arguments handling
            dots <- list(...)
            validDotsArgs(dots, names=c("maxIter","sandaft"))
            
            ## initializations ('ols' is always TRUE so not specified)
            rr <- aft <- uv <- FALSE
            aftfit <- rrfit <- uvfit <- "No fit was requested."

            if(data@has.censored) { type <- "aft" }
            else if(missing(type)) { type <- "rr" }

            type <- validFitType(type)

            maxIter <- if(is.null(dots$maxIter)) 100 else dots$maxIter 
            validNumeric(maxIter, positive=TRUE, integer=TRUE)

            if(type=="rr") { rr <- TRUE }
            else if(type=="aft") { aft <- TRUE }
            else if(type=="uv") { uv <- TRUE }
            ## End input arguments handling
            
            dfru <- data@dfru
            settings <- data@settings
            
            endptscale <- settings$endptscale
            grpnames <- settings$grpnames

            oldop <- options(contrasts=c("contr.treatment", "contr.poly"))
            on.exit(oldop, add=TRUE)

            validAft(type, dfru)
            
            ## Ordinary Least Squares is *always* fit
            olsfit <- if(endptscale=="log") {
              lm(log(endpt) ~ -1 + grpf, data=dfru)
            }
            else {
              lm(endpt ~ -1 + grpf, data=dfru)
            }
            olsfit$dfru <- dfru
            names(olsfit$coef) <- grpnames

            if(rr) {
              rrfit <-
                if(endptscale=="log") {
                  try(rlm(log(endpt) ~ -1 + grpf, data=dfru, method="MM",
                          maxit=maxIter, ...))
                }
                else {
                  try(rlm(endpt ~ -1 + grpf, data=dfru, method="MM",
                          maxit=maxIter, ...))
                }
              if(class(rrfit)[1]!="try-error") {
                ## that is, the number of iterations
                ## did not exceed specified limit maxIter
                ## and thus convergence occurred
                rrfit$dfru <- dfru
                names(rrfit$coef) <- grpnames

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
            else if(aft) {
              if(!is.null(sandaft <- dots$sandaft)) {
                validBoolean(sandaft)
              }
              else {
                sandaft <- TRUE
              }

              thesurvobject <- with(dfru,
                                    if(endptscale=="log") {
                                      survival::Surv(time=log(endpt1), time2=log(endpt2),
                                                     event=status, type="interval")
                                    }
                                    else {
                                      survival::Surv(time=endpt1, time2=endpt2,
                                                     event=status, type="interval")
                                    })
              aftfit <-  try(survreg(thesurvobject ~ -1 + grpf,
                                     data=dfru, dist="gaussian", maxiter=maxIter))

              
              if(class(aftfit)[1]!="try-error") {
                ## that is, the number of iterations
                ## did not exceed specified limit maxIter
                ## and thus convergence occurred
                aftfit$dfru <- dfru
                aftfit$Surv <- thesurvobject
                aftfit$maxIter <- maxIter 
                names(aftfit$coef) <- grpnames
                aftfit$sandaft <- sandaft
                if(sandaft) {
                  ## Manually compute the sandwich estimate since robust=TRUE
                  ## fails for survreg when cluster is omitted. The actual bug
                  ## appears with the application of rowsum, which is not
                  ## really needed since we assume all subjects are
                  ## independent (i.e. no clusters).
                  aftfit$naive.var <- aftfit$var
                  aftfit$var <- crossprod(resid(aftfit, "dfbeta"))
                }
              }
              else {
                warning(cgMessage("The Accelerated Failure Time (aft) fit did not",
                                  "converge in the specified number of",
                                  "iterations. You may want to try again with an",
                                  "increased value for the maxIter argument.",
                                  warning=TRUE))
                aftfit <- paste("The AFT fit did not converge in the specified",
                                "number of iterations.")
              }
            }
            else if(uv) {
              uvfit <- if(endptscale=="log") {
                gls(log(endpt) ~ -1 + grpf, data=dfru,
                    weights=varIdent(form = ~ 1 | grpf))
              }
              else {
                gls(endpt ~ -1 + grpf, data=dfru,
                    weights=varIdent(form = ~ 1 | grpf))
              }
              uvfit$dfru <- dfru
              names(uvfit$coef) <- grpnames
            }

            returnObj <- new("cgOneFactorFit",
                             olsfit=olsfit,
                             rrfit=rrfit,
                             aftfit=aftfit,
                             uvfit=uvfit,
                             settings=settings)

            ## return
            returnObj
          })


setMethod("print", "cgOneFactorFit",
          print.cgOneFactorFit <-
          function(x, title=NULL, endptname=NULL,...) {
            ##
            ## PURPOSE: Simple print version of cgOneFactorFit
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
            aftfit <- x@aftfit
            uvfit <- x@uvfit

            settings <- x@settings

            ols <- rr <- uv <- aft <- FALSE  ## Initializations
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
            if(class(olsfit)[1]=="lm" && (model!="rronly") && !aft && !uv) {
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
            if(endptname!="") { cat(paste("Endpoint:", endptname, "\n")) }

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

            if(aft) {
              cat("\nAccelerated Failure Time Model Fit\n")
              print(aftfit, ...)
            }

            if(uv) {
              cat("\nUnequal Variances Model Fit\n")
              print(uvfit, ...)
            }

            invisible()  
          }
          )


setMethod("show", "cgOneFactorFit", function(object) print(object))

setMethod("showObj", "cgOneFactorFit",
          showObj.cgOneFactorFit <- function(object) showDefault(object))

setMethod("summary", "cgOneFactorFit",
          summary.cgOneFactorFit <-
          function(object, title=NULL, endptname=NULL, ...) {
            ##
            ## PURPOSE: Simple summary of cgOneFactorFit
            ## object. Echoes summary methods for individual
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
            
            olsfit <- object@olsfit
            rrfit <- object@rrfit
            aftfit <- object@aftfit
            uvfit <- object@uvfit

            settings <- object@settings

            ols <- rr <- uv <- aft <- FALSE  ## Initializations
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
            if(class(olsfit)[1]=="lm" && (model!="rronly") && !aft && !uv) {
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
            if(endptname!="") { cat(paste("Endpoint:", endptname, "\n")) }

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

            if(aft) {
              cat("\nAccelerated Failure Time Model Fit Summary\n")
              print(summary(aftfit, ...))
            }

            if(uv) {
              cat("\nUnequal Variances Model Fit Summary\n")
              print(summary(uvfit, ...))
            }

            invisible()  

          }
          )


validAft <- function(type, dfru) {
  if(type=="aft" && ncol(dfru)!=5) {
    stop(cgMessage("An accelerated failure time (AFT) model",
                   "cannot be fit as requested (type=\"aft\")",
                   "since the data frame does not seem to have",
                   "a censored status column or the required format.",
                   seeHelpFile("CGOneFactorFit")))
  }
  return(TRUE)
}

validFitType <- function(type) {
  x <- try(match.arg(type, c("ols","rr","aft","uv")))
  if(class(x)=="try-error") {
    stop(cgMessage("The type argument needs to evaluate to one",
                   "of: \"ols\", \"rr\", \"aft\", \"uv\".",
                   seeHelpFile("CGOneFactorFit")))
  }
  else return(x)
}

