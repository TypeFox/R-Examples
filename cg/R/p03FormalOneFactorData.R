## $Id: p03FormalOneFactorData.R 6053 2015-02-22 20:23:45Z bpikouni $

## One-Factor Unpaired Groups Case

## Formal Analysis methods for One-Factor Unpaired Groups Data

setClass("cgOneFactorGlobalTest",
         representation(ols.gpval="numericOrNULL", rr.gpval="numericOrNULL",
                        aft.gpval="numericOrNULL", uv.gpval="numericOrNULL",
                        settings="list"),
         prototype(ols.gpval=NULL, rr.gpval=NULL, aft.gpval=NULL,
                   uv.gpval=NULL,
                   settings=list()))

setMethod("globalTest", "cgOneFactorFit",
          globalTest.cgOneFactorFit <- 
          function(fit, display="print", ...) {
            ##
            ## PURPOSE: Derive global F-test p-values
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names="model")
            
            ## initializations
            aft <- ols <- rr <- uv <- FALSE
            rr.gpval <- ols.gpval <- aft.gpval <- uv.gpval <- NULL
            ##
            rrfit <- fit@rrfit
            olsfit <- fit@olsfit
            aftfit <- fit@aftfit
            uvfit <- fit@uvfit

            settings <- fit@settings

            modelarg <- getDotsArgName(dots, "model")
            if(!is.na(modelarg)) {
              model <- eval(parse(text=paste("dots$", modelarg, sep="")))
              model <- validArgMatch(model, choices=c("both", "olsonly","rronly"))
            }
            else {
              model <- "both"
            }

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

            if(rr || ols) {
              df.residual <- olsfit$df.residual
              df.grpf <- length(settings$grpnames) - 1
              upform <- update(formula(olsfit), . ~ grpf)
            }

            dfru <- olsfit$dfru

            if(ols) {
              ## calculations mimic print.summary.lm code snippet
              olsfit <- lm(upform, data=dfru)
              ols.gpval <- with(summary(olsfit), pf(fstatistic[1], 
                                                    fstatistic[2],
                                                    fstatistic[3], lower.tail
                                                    = FALSE))
            }
            
            if(rr) {
              ## if only 2 groups, we rely on the comparisonstable method for computing
              ## the pvalue
              if(df.grpf==1) {
                rr.gpval <- pairwisecompsmatrix(comparisonsTable(fit,
                                                                 type="pairwisereflect")@rr.comprs,
                                                grpnames=settings$grpnames)[3, 3]
              }
              else {
                ## As of current writing,
                ## The MASS summary.rlm does not compute or facilitate construction of a
                ## global F, or equivalently, an R^2. For now we use a placeholder method
                ## based on ad-hoc (witchcraft) development, which is essentially to
                ## re-fit a linear model with lm() and weights from the RR fit.
                olsw.rrfit <- lm(upform, data=dfru, weights=rrfit$w)
                rr.Fobs <- { (1 / (1 / summary.lm(olsw.rrfit)$r.squared - 1))*
                               (df.residual/df.grpf) }
                rr.gpval <- 1 - pf(rr.Fobs, df.grpf, df.residual)
              }
            }
            else if(aft) {
              upform <- update(formula(aftfit), . ~ grpf)
              aftfit <- survreg(upform, data = dfru, 
                                dist = "gaussian", maxiter = aftfit$maxIter)
              ## calculations mimic print.survreg code snippet
              chisq <- 2 * diff(aftfit$loglik)
              df <- sum(aftfit$df) - aftfit$idf
              aft.gpval <- 1 - pchisq(chisq, df)
            }
            else if(uv) {
              upform <- update(formula(uvfit), . ~ grpf)
              uvfit <- gls(model = upform, data = dfru,
                           weights = varIdent(form = ~1 | grpf))
              uv.gpval <- anova(uvfit)$"p-value"[2]
            }

            x <- new("cgOneFactorGlobalTest",
                     ols.gpval=ols.gpval, rr.gpval=rr.gpval, 
                     aft.gpval=aft.gpval, uv.gpval=uv.gpval,
                     settings=settings)
            if(display=="print") {
              print(x)
            }
            else if (display=="show"){
              showDefault(x)
            }
            ## else show nothing
            invisible(x)
          })


setMethod("print", "cgOneFactorGlobalTest",
          print.cgOneFactorGlobalTest <-
          function(x, title=NULL, endptname=NULL, ...) {
            ##
            ## PURPOSE: Semi-formatted print version of Global test p-value
            ## Three digits are used, and anything smaller than or equal to
            ## 0.0005 gets "< 0.001"
            ## 
            ## NOTE: Had to use x as an argument because of the system defined
            ## generic. I would have preferred to use object; hence the first
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

            olsgrps <- object@ols.gpval
            rrgrps <- object@rr.gpval
            aftgrps <- object@aft.gpval
            uvgrps <- object@uv.gpval

            ols <- rr <- uv <- aft <- FALSE  ## Initializations
            if(!is.null(aftgrps)) {
              aft <- TRUE
              validArgModel(...)
            }
            else if(!is.null(uvgrps)) {
              uv <- TRUE
              validArgModel(...)              
            }
            if(!is.null(rrgrps) && model!="olsonly" && !aft && !uv) {
              rr <- TRUE
            }
            if(!is.null(olsgrps) && (model!="rronly") && !aft && !uv) {
              ols <- TRUE
            }

            settings <- object@settings
            alpha <- settings$alpha
            
            if(is.null(title)) {
              title <- paste("Global Test P-value of", settings$analysisname) 
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
              cat("\nClassical Least Squares Model Fit:", fmtPvalue(olsgrps), "\n")
            }
            
            if(rr) {
              cat("\nResistant & Robust Model Fit:", fmtPvalue(rrgrps), "\n")
            }

            if(aft) {
              cat("\nAccelerated Failure Time Model Fit:", fmtPvalue(aftgrps), "\n")
            }

            if(uv) {
              cat("\nUnequal Variances Model Fit:", fmtPvalue(uvgrps), "\n")
            }

            invisible()
            
          })

setMethod("show", "cgOneFactorGlobalTest",
          show.cgOneFactorGlobalTest <- function(object) showDefault(object))

## Comparisons Tables
comparisons <- function(estimates,
                        varcovmatrix, errordf=Inf,
                        endptscale, mcadjust=FALSE,
                        alpha=0.05, 
                        type="pairwisereflect",
                        contrastmatrix=NULL, n,
                        offset=NULL, cnames="derive",
                        analysisname="", endptname="",
                        digits=NULL,
                        addpct=FALSE,
                        display="print") {
  ## 
  ## PURPOSE: Function for computing comparisons results
  ## and placing in a listed form.
  ##
  ## Input argument checking
  if(!hasArg(estimates)) reportInvalidArg("estimates")
  if(!hasArg(varcovmatrix)) reportInvalidArg("varcovmatrix")
  validAlpha(alpha)
  validBoolean(mcadjust)
  validBoolean(addpct)
  validErrorDf(errordf, varcovmatrix, n, mcadjust)
  validComparisonType(type, errordf)
  validCnames(type, cnames, nrow(contrastmatrix))
  estimates <- validEstimates(estimates)
  endptscale <- validArgMatch(endptscale, c("log","original"))
  validAddPct(addpct, endptscale)

  display <- validArgMatch(display, c("print","none","show"))
  
  ngrps <- length(estimates)
  grpnames <- names(estimates)

  if(is.null(contrastmatrix)) {
    L <- contrastMatrix(list(grpnames), type=type)
    ## if type="allgroupstocontrol", then the first grpname is set to be
    ## the reference group level.
  }
  else if(type=="custom") {
    L <- contrastmatrix
  }

  ## If only two groups were analyzed, then multiplicity adjustment has
  ## no effect. So we ensure mcadjust is set to FALSE
  if(ngrps==2) { mcadjust <- FALSE }

  if(mcadjust) {
    set.seed(17)
    multcompInform()

    glhtobj <- list(object=NULL, linfct=L, coef=estimates,
                    vcov=varcovmatrix,
                    type="user-defined", alternative="two.sided",
                    df=errordf, 
                    rhs=rep(0, nrow(L)))
    
    class(glhtobj) <- "glht"
    mcpci.fit <- confint(glhtobj, level=1-alpha)
    mcptest.fit <- summary(glhtobj)

    estdiff <- mcpci.fit$confint[,"Estimate"]
    sediff <- mcptest.fit$test$sigma
    lowerci <- mcpci.fit$confint[,"lwr"]
    upperci <- mcpci.fit$confint[,"upr"]
    pval <-  mcptest.fit$test$pvalues
    
  }
  else if(!mcadjust) {
    estdiff <- L %*% estimates
    sediff <- sqrt(diag(L %*%  varcovmatrix %*% t(L)))
    if(errordf=="approx") {
      ## Apply Satterthwaite approximation for special case of pairwise
      ## or all groups to control comparisons (-1 vs 1 contrasts)
      incidence <- t(apply(L, 1, function(x) { x!=0 }))
      thestnderrs <- sqrt(diag(varcovmatrix))
      errordf <- apply(incidence, 1,
                       function(x, thestnderrs, n) {
                         sepair <- thestnderrs[x]
                         npair <- n[x]
                         sese <- sqrt(sepair[1]^2 + sepair[2]^2)
                         return(sese^4/(sepair[1]^4/(npair[1] - 1) +
                                        sepair[2]^4/(npair[2] - 1)))
                       }, thestnderrs=thestnderrs, n=n)
    }
    
    tcrit <- qt(1 - alpha/2, errordf)
    lowerci <- estdiff - tcrit * sediff
    upperci <- estdiff + tcrit * sediff
    pval <- 2 * (1 - pt(abs(estdiff/sediff), errordf))
  }
  
  mcp <- data.frame(estimate=estdiff, se=sediff,
                    lowerci=lowerci, upperci=upperci,
                    pval=pval)

  ## Log Case:
  ## Transformations to Percent Difference Expression Scale
  ## And we need to specially handle the asymmetry of percent change,
  ## since decreases cannot be greater than 100%
  ## also note that if offset is a valid number,
  ## some calculations will be replaced. 
  if(endptscale=="log") {
    ## Focus only on the percent differences for now
    logscalest <- mcp[,1]

    if(is.null(offset)) {
      mcp[, 1] <- 100 * ( exp(logscalest) - 1 ) # Point Estimate
      mcp[, 2] <- 100 * ( exp(logscalest) * mcp[, 2] ) # Std Err
      mcp[, 3:4] <- 100 * ( exp(mcp[, 3:4]) - 1) # Confidence Limits
    }
    else { ## offset
      validNumeric(offset, positive=TRUE)
      mA <- exp(mcp$meanA) - offset
      mB <- exp(mcp$meanB) - offset
      mAmBratio <- mA/mB
      correctionFactor <- (mA * (mB + offset))/(mB * (mA + offset))
      
      mcp[, 1] <- 100 * ( mAmBratio - 1 ) # Point Estimate
      mcp[, 2] <- 100 * ( mAmBratio * mcp[, 2] ) # Std Err
      mcp[, 3:4] <- 100 * ( correctionFactor*exp(mcp[, 3:4]) - 1) # Confidence Limits
    }
  }

  ## Fetch, arrange, and add individual component estimates
  if(type!="custom") {
    complist <- strsplit(row.names(L), split=" vs. ", fixed=TRUE)
    compA <- sapply(complist, function(x) x[1])
    compB <- sapply(complist, function(x) x[2])
    
    compAlist <- as.list(compA)
    compBlist <- as.list(compB)
    
    compAindx <- unlist(sapply(compAlist,
                               function(x) which(grpnames==x)))
    compBindx <- unlist(sapply(compBlist,
                               function(x) which(grpnames==x)))
    sevec <- sqrt(diag(varcovmatrix))
    
    mcp$meanA <- estimates[compAindx]
    mcp$seA <- sevec[compAindx]
    mcp$meanB <- estimates[compBindx]
    mcp$seB <- sevec[compBindx]
  }
  else if(type=="custom") {
    L.A <- t(apply(L, 1, function(x) ifelse(x > 0, x, 0)))
    ## For the B side of the contrast, we need to remove the minus sign
    L.B <- t(apply(L, 1, function(x) ifelse(x < 0, -x, 0)))
    mcp$meanA <- L.A %*% estimates
    mcp$seA <- sqrt(diag(L.A %*%  varcovmatrix %*% t(L.A)))
    mcp$meanB <- L.B %*% estimates
    mcp$seB <- sqrt(diag(L.B %*%  varcovmatrix %*% t(L.B)))
    if(length(cnames)==1 && cnames=="derive") {
      names.ests <- names(estimates)
      names.A <- apply(L.A, 1,
                       function(x) {
                         indx <- (x > 0)
                         coef.ests <- x[indx]
                         coef.ests <- ifelse(coef.ests%%1!=0,
                                             signif(coef.ests, 3), coef.ests) 
                         label <- paste(paste(coef.ests, names.ests[indx], sep="*"),
                                        collapse="+")
                         if(length(coef.ests) > 1) {
                           label <- paste("(", label, ")", sep="")
                         }
                         gsub("1\\*","", label)
                       })
      names.B <- apply(L.B, 1,
                       function(x) {
                         indx <- (x > 0)
                         coef.ests <- x[indx]
                         coef.ests <- ifelse(coef.ests%%1!=0,
                                             signif(coef.ests, 3), coef.ests) 
                         label <- paste(paste(coef.ests, names.ests[indx], sep="*"),
                                        collapse="+")
                         if(length(coef.ests) > 1) {
                           label <- paste("(", label, ")", sep="")
                         }
                         gsub("1\\*","", label)
                       })
      rownames(mcp) <- paste(names.A, names.B, sep=" vs. ")
    }
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
  
  if(cnames[1]!="derive") { row.names(mcp) <- cnames }

  if(endptscale=="original" && addpct==TRUE) {
    mcp$pctdiff <- 100*(mcp$meanA / mcp$meanB - 1)
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
      if(!is.null(digits)) {
        fmt.mcp$meanA <- fround(fmt.mcp$meanA, digits=digits)
        fmt.mcp$seA <- fround(fmt.mcp$seA, digits=digits)
        fmt.mcp$meanB <- fround(fmt.mcp$meanB, digits=digits)
        fmt.mcp$seB <- fround(fmt.mcp$seB, digits=digits)
      }
      { names(fmt.mcp)[is.element(names(fmt.mcp),
                                  c("meanA","seA",
                                    "meanB","seB"))] <- c("geomeanA","seA",
                                                          "geomeanB","seB") }
      { names(mcp)[is.element(names(mcp), c("meanA","seA",
                                            "meanB","seB"))] <-
                                              c("geomeanA","seA",
                                                "geomeanB","seB") }
    }
    else { ## Simple Differences
      pctdiff <- mcp$pctdiff
      if(!is.null(digits)) fmt.mcp <- fround(fmt.mcp, digits)
      fmt.mcp$pval <- pvalfmt
      if(endptscale=="original" && addpct==TRUE) {
        fmt.mcp$pctdiff <- fmtPercent(pctdiff)
      }
    }

    cat(paste(round(100*(1-alpha), 0), "% Confidence ",
              "(alpha of ", alpha, ")",
              if(mcadjust) {", Multiplicity Adjusted"},
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
  invisible(mcp)

}


setClass("cgOneFactorComparisonsTable",
         representation(ols.comprs="dataframeMatrixOrNULL",
                        rr.comprs="dataframeMatrixOrNULL",
                        aft.comprs="dataframeMatrixOrNULL",
                        uv.comprs="dataframeMatrixOrNULL",
                        settings="list"),
         prototype(ols.comprs=NULL, 
                   rr.comprs=NULL, 
                   aft.comprs=NULL, 
                   uv.comprs=NULL, 
                   settings=list()))

setMethod("comparisonsTable", "cgOneFactorFit",
          comparisonsTable.cgOneFactorFit <-
          function(fit,  
                   ## mcadjust=FALSE, 
                   type="pairwisereflect",
                   ## contrastmatrix=NULL,
                   ## refgrp=NULL,
                   alpha=0.05, addpct=FALSE,
                   display="print", ...) {
            ##
            ## PURPOSE: Produce multiple comparisons of groups based
            ## on a one-factor unpaired samples fit.
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names=c("mcadjust", "contrastmatrix", "refgrp","model"))
            type <- validArgMatch(type, c("pairwisereflect","pairwise",
                                          "allgroupstocontrol", "custom"))
            validAlpha(alpha)
            validBoolean(addpct)

            display <- validArgMatch(display, c("print","none","show"))

            mcadjustarg <- getDotsArgName(dots, "mcadjustarg")
            if(!is.na(mcadjustarg)) {
              mcadjust <- eval(parse(text=paste("dots$", mcadjustarg, sep="")))
              validBoolean(mcadjust)
            }
            else {
              mcadjust <- FALSE
            }
            
            contrastmatrixarg <- getDotsArgName(dots, "contrastmatrix")
            if(is.na(contrastmatrixarg)) {
              contrastmatrix <- 1
              contrastmatrix <- as.null(contrastmatrix)			
            }
            
            refgrparg <- getDotsArgName(dots, "refgrp")
            settings <- fit@settings
            if(!is.na(refgrparg)) {
              refgrp <- eval(parse(text=paste("dots$", refgrparg, sep="")))
              if(!is.null(refgrp)) { 
                refgrp <- validArgMatch(refgrp, settings$grpnames) 
              }
              if(type!="allgroupstocontrol") {
                stop(cgMessage("If the refgrp= argument is specified",
                   "then the type= argument needs to be set to",
                   "\"allgroupstocontrol\"."))
              }
            }
            else refgrp <- settings$refgrp
            
            modelarg <- getDotsArgName(dots, "model")
            if(!is.na(modelarg)) {
              model <- eval(parse(text=paste("dots$", modelarg, sep="")))
              model <- validArgMatch(model, choices=c("both", "olsonly","rronly"))
            }
            else {
              model <- "both"
            }
            
            ## initializations
            aft <- ols <- rr <- uv <- FALSE
            aft.comprs <- ols.comprs <- rr.comprs <- uv.comprs <- NULL
            ##
            ## settings <- fit@settings
            endptscale <- settings$endptscale
            rrfit <- fit@rrfit
            olsfit <- fit@olsfit
            aftfit <- fit@aftfit
            uvfit <- fit@uvfit
            grpnames <- settings$grpnames
            offset <- settings$addconstant
            
            validAddPct(addpct, endptscale)

            if(class(aftfit)[1]=="survreg") {
              aft <- TRUE
              validArgModel(...)
            }
            else if(class(uvfit)[1]=="gls") {
              uv <- TRUE
              validArgModel(...)
              thens <- with(uvfit$dfru, sapply(split(endpt, grpf), length))
            }

            if(class(rrfit)[1]=="rlm" && model!="olsonly" && !aft && !uv) {
              rr <- TRUE
            }
            if(class(olsfit)[1]=="lm" && model!="rronly" && !aft && !uv) {
              ols <- TRUE
              if(!rr) model <- "olsonly"
            }

            numberofgrps <- length(grpnames)
            grpnamesindex <- 1:numberofgrps
            if(!is.null(refgrp) && type=="allgroupstocontrol") {
              grpnamesindex <- c(which(grpnames==refgrp), which(grpnames!=refgrp))
            }
            
            if(rr || ols) {
              df.residual <- olsfit$df.residual
            }

            if(ols) {
              olsestimates <- olsfit$coef
              varcov.olsfit <- vcov(olsfit)
              ols.comprs <- comparisons(olsestimates[grpnamesindex],
                                        varcov.olsfit[grpnamesindex, grpnamesindex],
                                        errordf=df.residual,
                                        endptscale=endptscale,
                                        mcadjust=mcadjust,
                                        alpha=alpha,
                                        type=type,
                                        contrastmatrix=contrastmatrix,
                                        offset=offset, addpct=addpct, display="none")
              if(mcadjust) { multcompDone("Classical Least Squares") }
            }
            
            if(rr) {
              rrestimates <- rrfit$coef
              summ.rrfit <- summary(rrfit, method="XtWX", ...)
              stddev <- summ.rrfit$stddev
              varcov.rrfit <- (stddev^2) * summ.rrfit$cov.unscaled

              rr.comprs <- comparisons(rrestimates[grpnamesindex],
                                       varcov.rrfit[grpnamesindex, grpnamesindex],
                                       errordf=df.residual,
                                       endptscale=endptscale,
                                       mcadjust=mcadjust,
                                       alpha=alpha,
                                       type=type,
                                       contrastmatrix=contrastmatrix,
                                       offset=offset, 
                                       addpct=addpct, display="none")
              if(mcadjust) { multcompDone("Resistant & Robust") }
            }              
            
            if(aft) {
              aftestimates <- aftfit$coef
              ## estscaleindex <- length(aftestimates) + 1
              df.residual <- aftfit$df.residual
              varcov.aftfit <- vcov(aftfit)
              
              aft.comprs <- comparisons(aftestimates[grpnamesindex],
                                        varcov.aftfit[grpnamesindex, grpnamesindex],
                                        errordf=df.residual,
                                        endptscale=endptscale,
                                        mcadjust=mcadjust,
                                        alpha=alpha,
                                        type=type,
                                        contrastmatrix=contrastmatrix,
                                        offset=offset, display="none")
              if(mcadjust) { multcompDone("Accelerated Failure Time") }
            }

            if(uv) {
              uvestimates <- uvfit$coef
              varcov.uvfit <- vcov(uvfit)
              
              uv.comprs <- comparisons(uvestimates[grpnamesindex],
                                       varcov.uvfit[grpnamesindex, grpnamesindex],
                                       errordf="approx",
                                       endptscale=endptscale,
                                       mcadjust=mcadjust,
                                       alpha=alpha,
                                       type=type,
                                       contrastmatrix=contrastmatrix,
                                       n=with(uvfit$dfru,
                                         sapply(split(endpt, grpf),
                                                length)[grpnamesindex]),
                                       offset=offset, display="none")
            }

            ## slightly different from the original fit object settings slot
            settings <- list(endptscale=settings$endptscale,
                             analysisname=settings$analysisname,
                             endptname=settings$endptname,
                             endptlabel=makeEndptLabel(settings$endptname,
                               settings$endptunits), 
                             alpha=alpha,
                             mcadjust=mcadjust,
                             digits=settings$digits,
                             stamps=settings$stamps,
                             type=switch(type,
                               pairwise="Pairwise Comparisons (halfset)",
                               pairwisereflect="All Pairwise Comparisons",
                               allgroupstocontrol="All Groups versus Control",
                               custom=""),
                             sandaft=if(aft) TRUE else FALSE,
                             addpct=addpct)

            returnObj <- new("cgOneFactorComparisonsTable",
                             ols.comprs=ols.comprs,
                             rr.comprs=rr.comprs,
                             aft.comprs=aft.comprs,
                             uv.comprs=uv.comprs,
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

pairwisecompsmatrix <- function(comps, grpnames) {
  ## 
  ## PURPOSE: Function for rearranging pairwise comparisons results
  ## and placing in a matrix. The output matrix is designed to be
  ## used as input to formatting or further processing 
  ##
  ## comps is assumed to be a data frame such as
  ## created by the comparisons function
  grpnames <- validGrpNames(grpnames)
  ngrps <- length(grpnames)
  if(nrow(comps)!=ngrps*(ngrps-1)) {
    stop(cgMessage("The number of rows in the comps dataframe",
                   "needs to have entries for all possible",
                   "pairwise comparisons, i.e. the number of groups",
                   "times the number of groups - 1"))
  }
  ##ngrps <- round(uniroot(function(x) x^2 - x - nrow(comps),
  ##                       c(2, 1000))$root, 0)
  indx <- matrix(1:(ngrps^2) - as.vector(col(diag(ngrps))), ncol=ngrps)
  mcp.indx <- indx[lower.tri(indx)]
  mcp.opp.indx <- setdiff(1:nrow(comps), mcp.indx)

  mcp <- comps[mcp.indx,]
  mcp.opp <- comps[mcp.opp.indx,]

  ## The combined rows of mcp and mcp.opp
  ## pertain to all possible unique pairwise comparisons
  ## amongst the groups.  So we create a matrix for each statistical
  ## quantity that is of dimension ngrps by ngrps, putting the two
  ## triangular parts that come from mcp and mcp.opp together

  ## Construct list of each individual quantity
  mcp.matrices <- lapply(seq(along=mcp),
                         function(i, x, y, ngrps) {
                           z <- matrix(NA, nrow=ngrps, ncol=ngrps)
                           z[lower.tri(z)] <- x[[i]]
                           z[upper.tri(z)] <- y[[i]]
                           return(z)
                         }, x=mcp, y=mcp.opp, ngrps=ngrps)
  
  names(mcp.matrices) <- c("diff","se","lcl","ucl","pval")

  ## Now we have to combine the matrices into one big table
  ## which is probably intended for import into something like 
  ## EXCEL for formatting eventually.  This will involve
  ## interleaving rows and columns and adding some buffer sections to
  ## columns. A tricky part comes with the
  ## "lower CI and upper CI" part, because we want these to now be on the
  ## same row for a given comparison.  The dimension of the final table
  ## matrix will be 4*ngrps rows by 2*ngrps columns.
  mcpout.nrows <- 4*ngrps
  mcpout.ncols <- 2*ngrps
  mcp.out <- matrix(NA, nrow=mcpout.nrows, ncol=mcpout.ncols)
  ## fill-in the ucl's first
  mcp.out[rep( c(rep(FALSE, mcpout.nrows),
                 rep(c(rep(FALSE, 3), TRUE), ngrps)), ngrps ) ]  <- mcp.matrices$ucl

  ## second, the other quantities
  mcp.out[rep( c(rep(TRUE, mcpout.nrows),
                 rep(FALSE, mcpout.nrows) ),
              ngrps) ] <- do.call("rbind",
                                  lapply(mcp.matrices[c("diff","se",
                                                        "pval",
                                                        "lcl")],
                                         as.vector))
  ## lastly, matrix labels
  therownames <- character(length(grpnames)*4)
  thecolnames <- character(length(grpnames)*2)
  therownames[seq(1, length(therownames), by=4)] <- grpnames
  thecolnames[seq(1, length(thecolnames), by=2)] <- grpnames
  grpnameheaders <- list(therownames, thecolnames)
  dimnames(mcp.out) <- grpnameheaders
  
                                        # return
  mcp.out

}

setMethod("print", "cgOneFactorComparisonsTable",
          print.cgOneFactorComparisonsTable <-
          function(x, digits=NULL, title=NULL, endptname=NULL, ...) {
            ##
            ## PURPOSE: Semi-formatted print version of Comparisons Table
            ## 
            ## NOTE: Had to use x as an argument because of the system defined
            ## generic. I would have preferred to use object; hence the first
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

            olscomprs <- object@ols.comprs
            rrcomprs <- object@rr.comprs
            aftcomprs <- object@aft.comprs
            uvcomprs <- object@uv.comprs

            ols <- rr <- uv <- aft <- FALSE  ## Initializations
            if(!is.null(aftcomprs)) {
              aft <- TRUE
              validArgModel(...)
            }
            else if(!is.null(uvcomprs)) {
              uv <- TRUE
              validArgModel(...)              
            }
            if(!is.null(rrcomprs) && model!="olsonly" && !aft && !uv) {
              rr <- TRUE
            }
            if(!is.null(olscomprs) && (model!="rronly") && !aft && !uv) {
              ols <- TRUE
            }

            settings <- object@settings
            alpha <- settings$alpha
            mcadjust <- settings$mcadjust
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

            fmtdig <- function(x, diffmetric, digits) {
              pvalfmt <- fmtPvalue(x$pval)
              fmt.x <- x
              if(regexpr("Percent", diffmetric) > 0) {
                fmt.x$estimate <- fmtPercent(x$estimate)
                fmt.x$se <- fmtPercent(x$se)
                fmt.x$lowerci <- fmtPercent(x$lowerci)
                fmt.x$upperci <- fmtPercent(x$upperci)
                fmt.x$pval <- pvalfmt
                fmt.x$meanA <- fround(x$meanA, digits=digits)
                fmt.x$seA <- fround(x$seA, digits=digits)
                fmt.x$meanB <- fround(x$meanB, digits=digits)
                fmt.x$seB <- fround(x$seB, digits=digits)
                { names(fmt.x)[is.element(names(fmt.x), c("meanA","seA",
                                                          "meanB","seB"))] <-
                                                            c("geomeanA","seA",
                                                              "geomeanB","seB") }
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
                        if(mcadjust) {", Multiplicity Adjusted"},
                        "\n",
                        sep=""))
            }

            curwidth <- getOption("width")
            on.exit(options(width=curwidth), add=TRUE)
            if(curwidth < 500) { options(width=500) }
            
            if(ols) {
              cat("\nClassical Least Squares Model Fit\n")
              informConfidence()
              print(fmtdig(olscomprs, diffmetric, digits), quote=FALSE)
            }
            
            if(rr) {
              cat("\nResistant & Robust Model Fit\n")
              informConfidence()
              print(fmtdig(rrcomprs, diffmetric, digits), quote=FALSE)
            }

            if(aft) {
              afttitle <- "\nAccelerated Failure Time Model Fit\n"
              if(settings$sandaft) {
                cat(paste(afttitle,
                          "with Sandwich Variance-Covariance Estimate\n",
                          sep=""))
              }
              else {
                cat(afttitle)
              }
              informConfidence()
              print(fmtdig(aftcomprs, diffmetric, digits), quote=FALSE)
            }

            if(uv) {
              cat("\nUnequal Variances Model Fit\n")
              informConfidence()
              print(fmtdig(uvcomprs, diffmetric, digits), quote=FALSE)
            }
            
            invisible()
          })

setMethod("show", "cgOneFactorComparisonsTable",
          show.cgOneFactorComparisonsTable <- function(object) showDefault(object))


## Error Bar Graph
errorbargraph <- function(estimates, centralvar,
                          critpoint,
                          endptscale="log",
                          analysisname="",
                          endptname="",
                          alpha=0.05,
                          digits=NULL,
                          approxstamp=FALSE,
                          titlestamp=TRUE,
                          offset=NULL, 
                          ticklabels=NULL, ...) {
  ##
  ## PURPOSE: Construct an error bar graph based on pairwise
  ## multiple comparisons (method of Andrews, Sarner, and Snee, 1980)
  ##
  ## argument handling
  estimates <- validEstimates(estimates)
  endptscale <- validArgMatch(endptscale, c("log","original"))
  validAlpha(alpha)
  validBoolean(approxstamp)
  validBoolean(titlestamp)

  options(warn=-1)
  curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010)
  options(warn=0)
  on.exit(par(curpar))

  numberofgrps <- length(estimates)
  grpnames <- names(estimates)

  errorbarlength <- critpoint*sqrt(2*centralvar)/2
  bardata <-  list(lower=estimates - errorbarlength,
                   upper=estimates + errorbarlength)
  logscale <- ifelse(endptscale=="log", TRUE, FALSE)
  rangedata <- if(logscale) exp(unlist(bardata)) else unlist(bardata)
  
  if(logscale) {
    parmar <- par(mar=c(5, 4, 4, 3) + 0.1)
    curpar$mar <- parmar$mar
    
    errbar(1:numberofgrps, estimates/log(10),
           yplus=bardata$upper/log(10),
           yminus=bardata$lower/log(10),
           ylab=endptname, xlab="",
           pch=16,
           xlim=c(0.5, numberofgrps + 0.5), axes=FALSE)
    mtext("log-spaced", side=2, line=2.25, cex=0.7)
    tickmarks <- setupAxisTicks(rangedata,
                                logscale=logscale, digits=digits,
                                offset=offset)
    if(!is.null(ticklabels)) {
      tickmarks <- makeTickMarks(ticklabels, tickmarks,
                                 offset=offset)
    }
    
    axis(2, at=log10(tickmarks), labels=names(tickmarks),
         cex.axis=0.8, adj=1, las=1)
    
    log10endpt <- unlist(bardata)/log(10)
    axis(4, at=pretty(log10endpt), pretty(log10endpt), cex.axis=0.7,
         adj=0, las=1, tck=-0.0075)
    mtext(side=4, text=if(is.expression(endptname)) {
      catCharExpr("Log10 scale of", endptname)
    } else { paste("Log10 scale of", endptname, sep=" ") },
          line=2, cex=0.7, adj=0) 
  }
  else {
    errbar(1:numberofgrps, estimates,
           yplus=bardata$upper,
           yminus=bardata$lower,
           ylab=endptname, xlab="",
           pch=16,
           xlim=c(0.5, numberofgrps + 0.5), axes=FALSE)
    tickmarks <- setupAxisTicks(rangedata,
                                logscale=logscale, digits=digits)
    if(!is.null(ticklabels)) {
      tickmarks <- makeTickMarks(ticklabels, tickmarks,
                                 offset=offset)
    }

    axis(2, at=tickmarks, labels=names(tickmarks),
         cex.axis=0.8, adj=1, las=1)
  }

  ## Axes Customization
  grpnameticksettings <- setupGrpNameTicks(grpnames, 1:numberofgrps)
  plotGrpNameTicks((grpnames), settings=grpnameticksettings)
  
  minmaxTicks(if(logscale) exp(unlist(bardata)) else
              unlist(bardata),
              theaxis="y", logscale=logscale, digits=digits,
              offset=offset)

  ## Annotations
  if(titlestamp) {
    title(main=paste("Error Bar Graph\n",
            analysisname, sep=""), line=2, cex.main=1.1)
    errorBarGraphStamp(alphapercent=100*alpha)
  }
  if(approxstamp) errorBarGraphApproximateStamp()
  
  box()

  invisible() 
}

setMethod("errorBarGraph", "cgOneFactorFit",
          errorBarGraph.cgOneFactorFit <-
          function(fit, mcadjust=FALSE, 
                   alpha=0.05, cgtheme=TRUE, device="single", ...) {
            ##
            ## PURPOSE: Construct an error bar graph based on pairwise
            ## multiple comparisons (method of Andrews, Sarner, and Snee, 1980)
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names=c("model", "ticklabels"))
            
            if(class(fit@uvfit)[1]=="gls" || class(fit@aftfit)[1]=="survreg") {
              stop(cgMessage("There is no errorBarGraph method",
                             "defined for a fitted model that allowed",
                             "unequal variances or censored observations.",
                             "You may wish to look",
                             "at the errorbargraph function",
                             "to construct a graph if your standard errors",
                             "are not too different."))
            }

            validBoolean(cgtheme)
            validAlpha(alpha)
            validBoolean(mcadjust)
            
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
              model <- validArgMatch(model, choices=c("both", "olsonly", "rronly"))
            }
            else {
              model <- "both"
            }

            device <- validArgMatch(device, c("single", "multiple", "ask"))

            settings <- fit@settings
            offset <- settings$addconstant
            
            rrfit <- fit@rrfit
            olsfit <- fit@olsfit
            
            digits <- settings$digits
            analysisname <- settings$analysisname
            endptlabel <- makeEndptLabel(settings$endptname, settings$endptunits)            
            endptscale <- settings$endptscale
            stamps <- settings$stamps

            alphapercent <- round(100*alpha, 0)
            confquantile <- 1 - alpha/2

            grpnames <- settings$grpnames
            numberofgrps <- length(grpnames)
            if(numberofgrps==2) { mcadjust <- FALSE }
            L <- contrastMatrix(list(grpnames), "pairwise")

            ols <- rr <- FALSE  ## Initializations
            if(class(rrfit)[1]=="rlm" && model!="olsonly") {
              rr <- TRUE
            }
            if(class(olsfit)[1]=="lm" && model!="rronly") {
              ols <- TRUE
            }

            df.residual <- olsfit$df.residual
            thetitle <- "Error Bar Graph"
            dfru <- olsfit$dfru

            if(ols) {
              n <- with(dfru, tapply(endpt, grpf, length))
              ## Calculate the harmonic mean to account for unequal sample sizes
              nharmonic <- n ## default all n equal
              if(length(unique(n)) > 1) {
                nharmonic <- with(olsfit$dfru,
                                  length(n) / sum(tapply(endpt, grpf,
                                                         function(x) {
                                                           1/length(x)} )
                                                  ))
              }
              olsestimates <- olsfit$coef
              olscentralstderr <- summary(olsfit)$sigma/sqrt(nharmonic)
              olscentralvar <- olscentralstderr^2 

              olscritpoint <-
                if(!mcadjust) { qt(confquantile, df.residual)  }
                else {
                  set.seed(17)
                  multcompInform()

                  glhtobj <- list(object=NULL, linfct=L, coef=olsestimates,
                                  vcov=vcov(olsfit),
                                  type="user-defined", alternative="two.sided",
                                  df=df.residual, 
                                  rhs=rep(0, nrow(L)))
                  class(glhtobj) <- "glht"
                  ## assign this value:
                  attr(confint(glhtobj, level=1-alpha)$confint, "calpha")
                }
              if(mcadjust) { multcompDone("Classical Least Squares") }
            }

            if(rr) {
              rrestimates <- rrfit$coef
              summ.rrfit <- summary(rrfit, method="XtWX", correlation=FALSE)
              varcov.rrfit <- (summ.rrfit$stddev^2) * summ.rrfit$cov.unscaled

              ## a robust / resistant estimate of center of the std errors
              rrcentralvar <- rlm(diag(varcov.rrfit) ~ 1)[[1]]
              rrcentralstderr <- sqrt(rrcentralvar)

              rrcritpoint <-
                if(!mcadjust) { qt(confquantile, df.residual)  }
                else {
                  set.seed(17)
                  multcompInform()
                  
                  glhtobj <- list(object=NULL, linfct=L, coef=olsestimates,
                                  vcov=varcov.rrfit,
                                  type="user-defined", alternative="two.sided",
                                  df=df.residual, 
                                  rhs=rep(0, nrow(L)))
                  class(glhtobj) <- "glht"
                  ## assign this value:
                  attr(confint(glhtobj, level=1-alpha)$confint, "calpha")
                }
              if(mcadjust) { multcompDone("Resistant & Robust") }
              
            }

            if(rr && ols && is.element(model, "both") && device=="single") {
              ols.dfr <- data.frame(type=rep("Classical", numberofgrps),
                                    grpf=factorInSeq(grpnames), estimate=olsestimates,
                                    lower=olsestimates -
                                    olscritpoint*olscentralstderr/sqrt(2),
                                    upper=olsestimates +
                                    olscritpoint*olscentralstderr/sqrt(2))

              rr.dfr <- data.frame(type=rep("Resistant & Robust", numberofgrps),
                                   grpf=factorInSeq(grpnames), estimate=rrestimates,
                                   lower=rrestimates -
                                   rrcritpoint*rrcentralstderr/sqrt(2),
                                   upper=rrestimates +
                                   rrcritpoint*rrcentralstderr/sqrt(2))

              all.dfr <- rbind(ols.dfr, rr.dfr)
              thetitle <- paste(thetitle, "s", sep="")
              theapproxmsg <-
                if(n[1]!=nharmonic[1]) {
                  "Error bars are approximate."
                }
                else {
                  "Error bars in Resistant & Robust panel are approximate."
                }
              
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
                                   top=list(pad1=1, pad2=1),
                                   right=list(pad1=0.5, pad2=1)
                                   ))
              on.exit(trellis.par.set("axis.components", trellisparstg3),
                      add=TRUE)               

              ally <- unlist(all.dfr[, c("estimate","lower","upper")])

              if(endptscale=="log") {
                thegraph <- xyplot(Cbind(log10(exp(estimate)),
                                         log10(exp(lower)),
                                         log10(exp(upper))) ~ as.numeric(grpf)
                                   | type,
                                   data=all.dfr,
                                   digits=digits,
                                   panel = function(x, y, ...) {
                                     grpnameticksettings <- setupGrpNameTicks(grpnames,
                                                                              grplocation=
                                                                              1:length(grpnames),
                                                                              cexinit=0.8,
                                                                              cexthreshold=0.5,
                                                                              grid=TRUE)
                                     plotGrpNameTicks(grpnames, grpnameticksettings, grid=TRUE)
                                     if(is.element(panel.number(),
                                                   c(1, nlevels(all.dfr$type)))) {
                                       panelside <- ifelse(panel.number()==1,
                                                           "left", "right")
                                       ycex <- 0.6
                                       tickmarks <- setupAxisTicks(exp(ally),
                                                                   logscale=TRUE,
                                                                   grid=TRUE,
                                                                   digits=digits,
                                                                   offset=offset,
                                                                   ycex=ycex)
                                       if(!is.null(ticklabels)) {
                                         tickmarks <- makeTickMarks(ticklabels, tickmarks,
                                                                    offset=offset)
                                       }
                                       panel.axis(side=panelside,
                                                  at=log10(tickmarks),
                                                  labels=names(tickmarks),
                                                  tck=0.20, text.cex=ycex,
                                                  rot=0,
                                                  outside=TRUE)
                                     }
                                     panel.xYplot(x, y, subscripts=FALSE,
                                                  method="bars",
                                                  label.curves=FALSE, ...)
                                     yrange <- range(unwind(attr(y, "other")))
                                     panel.text(x=rep(0, 2), y=yrange,
                                                labels=paste("",
                                                  c(fround(min(10^yrange), digits),
                                                    fround(max(10^yrange), digits))),
                                                col="blue", adj=0, cex=0.6)
                                   },
                                   layout=c(2,1), aspect=1, pch=16,
                                   col="black",
                                   as.table=TRUE,
                                   xlim=c(0, numberofgrps + 0.5),
                                   ylim=rangeExtend(log10(exp(ally))),
                                   ylab=list(cex=0.7,
                                     label=endptlabel),
                                   xlab="",
                                   scales=list(x=list(at=1:numberofgrps,
                                                 labels=NULL, tck=c(0.15, 0), axs="i"),
                                     y=list(labels=NULL, tck=0)),
                                   main=list(label=paste(thetitle, "\n",
                                               settings$analysisname, sep=""), cex=1.1),
                                   par.strip.text=list(cex=0.7)
                                   )
              }
              else {
                thegraph <- xyplot(Cbind(estimate,
                                         lower,
                                         upper) ~ as.numeric(grpf) | type,
                                   data=all.dfr,
                                   digits=digits,
                                   panel = function(x, y, ...) {
                                     grpnameticksettings <- setupGrpNameTicks(grpnames,
                                                                              grplocation=
                                                                              1:length(grpnames),
                                                                              cexinit=0.8,
                                                                              cexthreshold=0.5,
                                                                              grid=TRUE)
                                     plotGrpNameTicks(grpnames, grpnameticksettings, grid=TRUE)
                                     if(is.element(panel.number(),
                                                   c(1, nlevels(all.dfr$type)))) {
                                       panelside <- ifelse(panel.number()==1,
                                                           "left", "right")
                                       ycex <- 0.6
                                       tickmarks <- setupAxisTicks(ally,
                                                                   logscale=FALSE,
                                                                   grid=TRUE,
                                                                   digits=digits,
                                                                   offset=offset,
                                                                   ycex=ycex)
                                       if(!is.null(ticklabels)) {
                                         tickmarks <- makeTickMarks(ticklabels, tickmarks,
                                                                    offset=offset)
                                       }

                                       panel.axis(side=panelside,
                                                  at=tickmarks,
                                                  labels=names(tickmarks),
                                                  tck=0.20, text.cex=ycex,
                                                  rot=0,
                                                  outside=TRUE)
                                     }
                                     panel.xYplot(x, y, subscripts=FALSE,
                                                  method="bars",
                                                  label.curves=FALSE, ...)
                                     yrange <- range(unwind(attr(y, "other")))
                                     panel.text(x=rep(0, 2), y=yrange,
                                                labels=paste("",
                                                  c(fround(min(yrange), digits),
                                                    fround(max(yrange), digits))),
                                                col="blue", adj=0, cex=0.6)
                                   },
                                   layout=c(2,1), aspect=1, pch=16,
                                   col="black",
                                   as.table=TRUE,
                                   xlim=c(0, numberofgrps + 0.5),
                                   ylim=rangeExtend(ally),
                                   ylab=list(cex=0.7,
                                     label=endptlabel),
                                   xlab="",
                                   scales=list(x=list(at=1:numberofgrps,
                                                 labels=NULL, tck=c(0.15, 0), axs="i"),
                                     y=list(labels=NULL, tck=0)),
                                   main=list(label=paste(thetitle, "\n",
                                               settings$analysisname, sep=""), cex=1.1),
                                   par.strip.text=list(cex=0.7)
                                   )
              }

              print(thegraph) ## , position=c(0,0,0.95,1))

              if(stamps) graphStampCG()
              errorBarGraphStamp(mcadjust, alphapercent, grid=TRUE)
              errorBarGraphApproximateStamp(grid=TRUE, msg=theapproxmsg)
              
            }

            else if((model=="olsonly" || (ols && !rr && model=="both")) &&
                    device=="single") {
              errorbargraph(olsestimates, olscentralvar,
                            olscritpoint,
                            endptscale,
                            analysisname,
                            endptname=endptlabel,
                            alpha=alpha,
                            digits=digits,
                            approxstamp={nharmonic[1]!=n[1]},
                            titlestamp=FALSE, offset=offset,
                            ticklabels=ticklabels)              
              
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Error Bar Graph, Classical analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              errorBarGraphStamp(mcadjust, alphapercent)

            }

            else if(model=="rronly" && rr && device=="single") {
              errorbargraph(rrestimates, rrcentralvar,
                            rrcritpoint,
                            endptscale,
                            analysisname,
                            endptname=endptlabel,
                            alpha=alpha,
                            digits=digits,
                            approxstamp=TRUE,
                            titlestamp=FALSE,
                            offset=offset,
                            ticklabels=ticklabels)
              
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Error Bar Graph, Resistant & Robust analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              errorBarGraphStamp(mcadjust, alphapercent)

            }

            else if(rr && ols &&
                    is.element(device, c("ask","multiple")) &&
                    is.element(model, "both")) { 
              
              device <- validArgMatch(device, c("multiple", "ask"))
              if(device=="ask") {
                op <- par(ask = TRUE)
                on.exit(par(op), add=TRUE)
              }

              errorbargraph(olsestimates, olscentralvar,
                            olscritpoint,
                            endptscale,
                            analysisname,
                            endptname=endptlabel,
                            alpha=alpha,
                            digits=digits,
                            approxstamp={nharmonic[1]!=n[1]},
                            titlestamp=FALSE, offset=offset,
                            ticklabels=ticklabels)
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Error Bar Graph, Classical analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              errorBarGraphStamp(mcadjust, alphapercent)

              if(device=="multiple") {
                if(!is.null(dots$model)) dots$model <- NULL ## since we only
                ## want dots arguments for trellis.device in next call
                do.call("cgDevice", c(list(new=TRUE), dots))
                cat(cgMessage("A new graphics device has been generated",
                              "to hold the Resistant & Robust",
                              "errorBarGraph version.",
                              "The Classical Least Squares version is on the previous",
                              "device.\n",
                              warning=TRUE))
              }
              
              errorbargraph(rrestimates, rrcentralvar,
                            rrcritpoint,
                            endptscale,
                            analysisname,
                            endptname=endptlabel,
                            alpha=alpha,
                            digits=digits,
                            approxstamp=TRUE,
                            titlestamp=FALSE,
                            offset=offset,
                            ticklabels=ticklabels)

              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Error Bar Graph, Resistant & Robust analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              errorBarGraphStamp(mcadjust, alphapercent)
              
            }

            else {
              stop(cgMessage("The chosen device and model arguments",
                             "are not compatible either with each other",
                             "or with the fitted model(s) in the fit object.",
                             seeHelpFile("errorBarGraph")))
            }
            
            invisible()
          }
          )

grpsummary <- function(estimates,
                       varcovmatrix,
                       n, errordf=Inf,
                       endptscale="log", mcadjust=FALSE,
                       alpha=0.05, 
                       offset=NULL, 
                       analysisname="", endptname="",
                       digits=NULL, display="print", ...) {
  ##
  ## PURPOSE: Function for computing group summary results.
  ## The output matrix is designed to be
  ## used as input to formatting or further processing
  validBoolean(mcadjust)
  n <- validN(n, estimates)
  estimates <- validEstimates(estimates)
  endptscale <- validArgMatch(endptscale, c("log","original"))
  
  validAlpha(alpha)
  display <- validArgMatch(display, c("print","none","show"))
  
  if(!is.numeric(errordf) || any(errordf < 1)) {
    stop(cgMessage("The errordf argument needs to be numeric",
                   "and (each value) set to be 1 or greater."))
  }

  ngrps <- length(estimates)
  grpnames <- names(estimates)

  ## If only two groups were analyzed, then multiplicity adjustment has
  ## no effect. So we ensure mcadjust is set to FALSE
  if(ngrps==2) { mcadjust <- FALSE }

  if(mcadjust) {
    if(length(errordf) > 1) {
      stop(cgMessage("The errordf can only be a single value",
                     "when mcadjust=FALSE.",
                     seeHelpFile("grpsummary")))
    }
    set.seed(17)
    multcompInform()

    glhtobj <- list(object=NULL, linfct=diag(ngrps),
                    coef=estimates,
                    vcov=varcovmatrix,
                    type="user-defined", alternative="two.sided",
                    df=errordf, 
                    rhs=rep(0, nrow(diag(ngrps))))
    class(glhtobj) <- "glht"
    mcpci.fit <- confint(glhtobj, level=1-alpha)
    mcptest.fit <- summary(glhtobj)

    se <- as.vector(mcptest.fit)$test$sigma
    lowerci <- mcpci.fit$confint[,"lwr"]
    upperci <- mcpci.fit$confint[,"upr"]
  }
  else if(!mcadjust) {
    se <- sqrt(diag(varcovmatrix))
    tcrit <- qt(1 - alpha/2, errordf)
    lowerci <- estimates - tcrit * se
    upperci <- estimates + tcrit * se
  }

  thetable <- data.frame(n=n,
                         estimate=estimates, se=se,
                         lowerci=lowerci, upperci=upperci)
  row.names(thetable) <- grpnames

  ## Handle log case if needed
  ## Transformations to Percent Difference Expression Scale
  if(endptscale=="log") {
    logscale.est <- thetable[, 2]
    logscale.sd <- thetable[, 3]
    
    thetable[, 2] <- exp(logscale.est) # Point Estimate
    thetable[, 3] <- exp(logscale.est) * logscale.sd  # Std Err
    if(any(logscale.sd > 0.50)) {
      warning(cgMessage("There is at least one",
                        "group standard error in the log scale",
                        "that exceeds 0.50, so the",
                        "estimated geometric mean standard",
                        "errors may be nonsensical and should be",
                        "cautiously regarded.")
              ) 
    }
    thetable[, 4:5] <- exp(thetable[, 4:5])  # Confidence Limits
  }

  if(!is.null(offset)) {
    validNumeric(offset, positive=TRUE)
    thetable[, c(2,4,5)] <- (thetable[, c(2,4,5)] - offset)
  }

  if(display=="print") {
    fmt.table <- thetable
    cat(paste("Group Summary Table",
              if(analysisname!="") paste("for", analysisname), "\n"))
    if(endptname!="") { cat(paste("Endpoint:", endptname, "\n")) }
    metric <- "Means"
    if(endptscale=="log") {
      metric <- paste("Geometric", metric)
    }
    else {
      metric <- paste("Arithmetic", metric)
    }
    cat(metric, "\n")

    if(!is.null(digits)) {
      fmt.table$estimate <- fround(fmt.table$estimate, digits)
      fmt.table$se <- fround(fmt.table$se, digits)
      fmt.table$lowerci <- fround(fmt.table$lowerci, digits)
      fmt.table$upperci <- fround(fmt.table$upperci, digits)
    }
    
    fmt.table$n <- fround(fmt.table$n, 0)

    cat(paste(round(100*(1-alpha), 0), "% Confidence ",
              "(alpha of ", alpha, ")",
              if(mcadjust) {", Multiplicity Adjusted"},
              "\n",
              sep=""))

    curwidth <- getOption("width")
    on.exit(options(width=curwidth), add=TRUE)
    if(curwidth < 500) { options(width=500) }
    
    print(fmt.table, quote=FALSE)
    
    return(invisible(thetable))
  }
  else if (display=="show") {
    showDefault(thetable)
  }

  ## else show nothing
  invisible(thetable)
}

setClass("cgOneFactorGrpSummaryTable",
         representation(ols.grps="dataframeOrNULL",
                        rr.grps="dataframeOrNULL",
                        aft.grps="dataframeOrNULL",
                        uv.grps="dataframeOrNULL",
                        settings="list"),
         prototype(ols.grps=NULL, rr.grps=NULL, aft.grps=NULL, uv.grps=NULL,
                   settings=list()))

setMethod("grpSummaryTable", "cgOneFactorFit",
          function(fit, 
                   mcadjust=FALSE,
                   alpha=0.05, display="print", ...) {
            ##
            ## PURPOSE: Summarize each group distribution based
            ## on the model fit.
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names="model")
            settings <- fit@settings
            
            validAlpha(alpha)
            validBoolean(mcadjust)
            display <- validArgMatch(display, c("print","none","show"))

            modelarg <- getDotsArgName(dots, "model")
            if(!is.na(modelarg)) {
              model <- eval(parse(text=paste("dots$", modelarg, sep="")))
              model <- validArgMatch(model, choices=c("both", "olsonly", "rronly"))
            }
            else {
              model <- "both"
            }
            ##
            
            ## initializations
            aft <- ols <- rr <- uv <- FALSE
            aft.grps <- ols.grps <- rr.grps <- uv.grps <- NULL
            ##
            endptscale <- settings$endptscale
            rrfit <- fit@rrfit
            olsfit <- fit@olsfit
            aftfit <- fit@aftfit
            uvfit <- fit@uvfit
            grpnames <- settings$grpnames
            offset <- settings$addconstant
            
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

            n <- with(olsfit$dfru, sapply(split(endpt, grpf), length))

            if(rr || ols) {
              df.residual <- olsfit$df.residual
            }

            if(ols) {
              olsestimates <- olsfit$coef
              varcov.olsfit <- vcov(olsfit)
              
              ols.grps <- grpsummary(olsestimates,
                                     varcov.olsfit,
                                     n=n, errordf=df.residual,
                                     endptscale=endptscale,
                                     mcadjust=mcadjust,
                                     alpha=alpha,
                                     offset=offset, display="none", ...)
              if(mcadjust) { multcompDone("Classical Least Squares") }
            }
            
            if(rr) {
              rrestimates <- rrfit$coef
              summ.rrfit <- summary(rrfit, method="XtWX")
              stddev <- summ.rrfit$stddev
              varcov.rrfit <- (stddev^2) * summ.rrfit$cov.unscaled

              rr.grps <- grpsummary(rrestimates,
                                    varcov.rrfit,
                                    n=n, errordf=df.residual,
                                    endptscale=endptscale,
                                    mcadjust=mcadjust,
                                    alpha=alpha,
                                    offset=offset, display="none", ...)
              if(mcadjust) { multcompDone("Resistant & Robust") }
            }

            if(aft) {
              aftestimates <- aftfit$coef
              df.residual <- aftfit$df.residual
              varcov.aftfit <- chop.matrix(vcov(aftfit))

              aft.grps <- grpsummary(aftestimates,
                                     varcov.aftfit,
                                     n=n, errordf=df.residual,
                                     endptscale=endptscale,
                                     mcadjust=mcadjust,
                                     alpha=alpha,
                                     offset=offset, display="none",...)
            }

            if(uv) {
              if(mcadjust) {
                stop(cgMessage("The mcadjust argument cannot be set",
                               "to TRUE when the fitted model allowed",
                               "for unequal variances."))
              }
              uv.grps <- grpsummary(uvfit$coef,
                                    vcov(uvfit),
                                    n=n, errordf=n-1,
                                    endptscale=endptscale,
                                    mcadjust=mcadjust,
                                    alpha=alpha,
                                    offset=offset, display="none",...)
            }

            settings <- list(endptscale=settings$endptscale,
                             analysisname=settings$analysisname,
                             endptname=settings$endptname,
                             alpha=alpha,
                             mcadjust=mcadjust,
                             digits=settings$digits,
                             sandaft=if(aft) TRUE else FALSE)

            returnObj <- new("cgOneFactorGrpSummaryTable",
                             ols.grps=ols.grps, rr.grps=rr.grps, 
                             aft.grps=aft.grps, uv.grps=uv.grps,
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


setMethod("print", "cgOneFactorGrpSummaryTable",
          print.cgOneFactorGrpSummaryTable <-
          function(x, digits=NULL, title=NULL, endptname=NULL, ...) {
            ##
            ## PURPOSE: Semi-formatted print version of Comparisons Table
            ## 
            ## NOTE: Had to use x as an argument because of the system defined
            ## generic. I would have preferred to use object; hence the first
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

            olsgrps <- object@ols.grps
            rrgrps <- object@rr.grps
            aftgrps <- object@aft.grps
            uvgrps <- object@uv.grps

            ols <- rr <- uv <- aft <- FALSE  ## Initializations
            if(!is.null(aftgrps)) {
              aft <- TRUE
              validArgModel(...)
            }
            else if(!is.null(uvgrps)) {
              uv <- TRUE
              validArgModel(...)              
            }
            if(!is.null(rrgrps) && model!="olsonly" && !aft && !uv) {
              rr <- TRUE
            }
            if(!is.null(olsgrps) && (model!="rronly") && !aft && !uv) {
              ols <- TRUE
            }

            settings <- object@settings
            alpha <- settings$alpha
            mcadjust <- settings$mcadjust
            
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
              title <- paste("Group Summary Table of", settings$analysisname) 
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

            metric <- "Means"
            if(settings$endptscale=="log") {
              metric <- paste("Geometric", metric)
            }
            else {
              metric <- paste("Arithmetic", metric)
            }
            cat(metric, "\n")

            fmtdig <- function(x, digits) {
              x$estimate <- fround(x$estimate, digits)
              x$se <- fround(x$se, digits)
              x$lowerci <- fround(x$lowerci, digits)
              x$upperci <- fround(x$upperci, digits)
              x$n <- fround(x$n, 0)
              return(x)
            }

            informConfidence <- function() {
              cat(paste(round(100*(1-alpha), 0), "% Confidence ",
                        "(alpha of ", alpha, ")",
                        if(mcadjust) {", Multiplicity Adjusted"},
                        "\n",
                        sep=""))
            }

            curwidth <- getOption("width")
            on.exit(options(width=curwidth), add=TRUE)
            if(curwidth < 500) { options(width=500) }
            
            if(ols) {
              cat("\nClassical Least Squares Model Fit\n")
              informConfidence()
              print(fmtdig(olsgrps, digits), quote=FALSE)
            }
            
            if(rr) {
              cat("\nResistant & Robust Model Fit\n")
              informConfidence()
              print(fmtdig(rrgrps, digits), quote=FALSE)
            }

            if(aft) {
              afttitle <- "\nAccelerated Failure Time Model Fit\n"
              if(settings$sandaft) {
                cat(paste(afttitle,
                          "with Sandwich Variance-Covariance Estimate\n",
                          sep=""))
              }
              else {
                cat(afttitle)
              }
              informConfidence()
              print(fmtdig(aftgrps,  digits), quote=FALSE)
            }

            if(uv) {
              cat("\nUnequal Variances Model Fit\n")
              informConfidence()
              print(fmtdig(uvgrps,  digits), quote=FALSE)
            }
            
            invisible()
          })

setMethod("show", "cgOneFactorGrpSummaryTable",
          show.cgOneFactorGrpSummaryTable <- function(object) showDefault(object))


comparisonsgraph <- function(compstable,
                             difftype,
                             analysisname="",
                             endptname="",
                             alpha=0.05,
                             digits=NULL,
                             titlestamp=TRUE,
                             explanation=TRUE,
                             wraplength=20,
                             cex.comps=0.7,
                             ticklabels=NULL,
                             ...) {
  ##
  ## PURPOSE: Plot Confidence Intervals of Comparisons;
  ## use horizontal error bars
  ##
  ## compstable needs to be a dataframe of proper format
  difftype <- validArgMatch(difftype, c("percent","amount","simple"))
  validAlpha(alpha)
  validBoolean(titlestamp)
  validBoolean(explanation)
  digits <- validArgDigits(digits)
  
  options(warn=-1)
  curpar <- par(new=FALSE, mgp=c(3,0.25,0), tck=-0.010, mar=c(5,7,4,2)+0.1)
  options(warn=0)
  on.exit(par(curpar), add=TRUE)

  numberofcomprs <- nrow(compstable)
  comprnames <- row.names(compstable)

  est <- compstable[, "estimate"]
  lower <- compstable[, "lowerci"]
  upper <- compstable[, "upperci"]
  compys <-  numberofcomprs:1

  if(difftype=="percent") {
    est <- log(pctToRatio(est))
    lower <- log(pctToRatio(lower))
    upper <- log(pctToRatio(upper))
  }
  
  allx <- c(lower, upper, est)

  ## We include zero because we always want the
  ## "no-difference" reference line to be visible
  xmin <- min(c(0, allx))
  xmax <- max(c(0, allx))
  
  if(difftype=="percent") {
    ## snippet from Hmisc's errbar function
    plot(est/log(10), compys, ylim = c(0.5, numberofcomprs + 0.5),
         xlim = c(xmin, xmax)/log(10),
         ylab="",
         xlab = "", 
         pch = 16, axes=FALSE, ...)
    mtext(side=1, text="log-spaced", line=4, cex=0.7)
    ycoord <- par()$usr[3:4]
    segments(lower/log(10), compys, upper/log(10), compys)
    smidge <- 0.015 * (ycoord[2] - ycoord[1])/2
    segments(lower/log(10), compys - smidge, lower/log(10), compys + smidge)
    segments(upper/log(10), compys - smidge, upper/log(10), compys + smidge)

    ## zero reference line
    abline(v=0, lty=3)
    
    ## Axes Customization
    log10endpt <- log10(exp(allx))
    xratioticks <- setupAxisTicks(exp(allx), ratio=TRUE, percent=TRUE,
                                  axis="x", digits=digits)

    if(!is.null(ticklabels)) {
      xratioticks <- makeTickMarks(ticklabels, xratioticks,
                                   percent=TRUE)
    }
    
    axis(1,at=log10(xratioticks),
         labels=names(xratioticks),
         tck=-.010, cex.axis=0.8)

    xlabchar <- "Percent Difference"
  }

  else { ## difftype = "amount" or "simple" (Simple Difference)
    plot(est, compys, ylim = c(0.5, numberofcomprs + 0.5),
         xlim = c(xmin, xmax),
         ylab="",
         xlab="",## paste("Difference in\n ", endptname, sep=""),
         pch=16, axes=FALSE, ...)
    ycoord <- par()$usr[3:4]
    segments(lower, compys, upper, compys)
    smidge <- 0.015 * (ycoord[2] - ycoord[1])/2
    segments(lower, compys - smidge, lower, compys + smidge)
    segments(upper, compys - smidge, upper, compys + smidge)

    ## zero reference line
    abline(v=0, lty=3)
    
    ## Axes Customization
    xdiffticks <- setupAxisTicks((allx), ratio=FALSE, difference=TRUE,
                                 logscale=FALSE,
                                 axis="x", digits=digits)
    if(!is.null(ticklabels)) {
      xdiffticks <- makeTickMarks(ticklabels, xdiffticks,
                                  percent=FALSE)
    }
    
    axis(1, at=xdiffticks,
         labels=names(xdiffticks),
         tck=-.010, cex=0.8)
    xlabchar <- "Difference"
  }

  mtext(side=1, text=xlabchar, line=2)
  if(is.expression(endptname)) {
    mtext(side=1, text=catCharExpr("in ", endptname), line=3)
  }
  else {
    mtext(side=1, text=paste("in", endptname), line=3)
  }
  
  ## Wrap long comparison labels
  for(i in seq(along=comprnames)) {
    ##if(nchar(comprnames[i]) > wraplength) {
    ##  comprnames[i] <- gsub("vs.","vs.\n", comprnames[i])
    ##}
    if(nchar(comprnames[i]) > wraplength) {
      comprnames[i] <- paste(strwrap(comprnames[i], width=wraplength,
                                     exdent=1), collapse="\n")
    }
  }
  axis(2, at=compys, labels=comprnames,
       adj=1, cex.axis=cex.comps, las=1)
  box()
  
  ## Titles
  par(curpar)
  if(explanation) {
    comparisonsGraphStamp(alphapercent=100*alpha)
  }
  if(titlestamp) {
    title(main=paste("Comparisons Graph\n",
            analysisname, sep=""), line=2, cex.main=1.1)    
  }

  ## min-max annotation
  bardata <- c(est, lower, upper)
  minmaxTicks(if(difftype=="percent") exp(bardata) else bardata,
              theaxis="x", percent=if(difftype=="percent") TRUE else FALSE,
              ratio=if(difftype=="percent") TRUE else FALSE,
              logscale=if(difftype=="percent") TRUE else FALSE,
              digits=digits)
  ## Add 0 
  text(y=par("usr")[3], x=0,
       labels="0\n", col="blue", adj=0, cex=0.7)

  invisible()

}

setMethod("comparisonsGraph", "cgOneFactorComparisonsTable",
          comparisonsGraph.cgOneFactorComparisonsTable <-
          function(compstable, cgtheme=TRUE, device="single",
                   wraplength=20, cex.comps=0.7, ...) {
            ##
            ## PURPOSE: create a graph of CI's on differences
            ## to judge comparisons
            ##
            ## Input arguments check
            dots <- list(...)
            validDotsArgs(dots, names=c("model", "ticklabels"))
            validBoolean(cgtheme)

            settings <- compstable@settings
            difftype <- if(settings$endptscale=="log") "percent" else "simple"
            alpha <- settings$alpha
            alphapercent <- round(100*alpha, 0)
            mcadjust <- settings$mcadjust
            analysisname <- settings$analysisname
            endptlabel <- settings$endptlabel
            stamps <- settings$stamps
            digits <- settings$digits

            rrcomprs <- compstable@rr.comprs
            olscomprs <- compstable@ols.comprs
            aftcomprs <- compstable@aft.comprs
            uvcomprs <- compstable@uv.comprs

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

            aft <- ols <- rr <- uv <- FALSE  ## Initializations
            if(!is.null(aftcomprs)) {
              aft <- TRUE
              validArgModel(...)
            }
            else if(!is.null(uvcomprs)) {
              uv <- TRUE
              validArgModel(...)              
            }
            if(!is.null(rrcomprs) && model!="olsonly" && !aft && !uv) {
              rr <- TRUE
            }
            if(!is.null(olscomprs) && (model!="rronly") && !aft && !uv) {
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

                                     ## ycoord <- par()$usr[3:4]
                                     ## smidge <- 0.015 * (ycoord[2] - ycoord[1])
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

              comparisonsGraphStamp(mcadjust, alphapercent, grid=TRUE,
                                    desc=settings$type)
              if(stamps) graphStampCG()
              
            }
            else if((model=="olsonly" || (ols && !rr && model=="both")) &&
                    device=="single") {
              comparisonsgraph(olscomprs, difftype, analysisname,
                               endptlabel, alpha, digits, titlestamp=FALSE,
                               explanation=FALSE,
                               wraplength=wraplength,
                               cex.comps=cex.comps,
                               ticklabels=ticklabels)
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Comparisons Graph, Classical analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              comparisonsGraphStamp(mcadjust, alphapercent,
                                    desc=settings$type)
            }
            
            else if((model=="rronly" || (!ols && rr && model=="both")) &&
                    device=="single") {
              ##            else if(model=="rronly" & rr & device=="single") {
              comparisonsgraph(rrcomprs, difftype, analysisname,
                               endptlabel, alpha, digits, titlestamp=FALSE,
                               explanation=FALSE,
                               wraplength=wraplength,
                               cex.comps=cex.comps,
                               ticklabels=ticklabels)
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Comparisons Graph, Resistant & Robust analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              comparisonsGraphStamp(mcadjust, alphapercent,
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
                               endptlabel, alpha, digits, titlestamp=FALSE,
                               explanation=FALSE,
                               wraplength=wraplength,
                               cex.comps=cex.comps,
                               ticklabels=ticklabels)
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Comparisons Graph, Classical analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              comparisonsGraphStamp(mcadjust, alphapercent,
                                    desc=settings$type)
              
              if(device=="multiple") {
                ## dots$model <- NULL ## avoid unused argument error
                
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
                               endptlabel, alpha, digits, titlestamp=FALSE,
                               explanation=FALSE,
                               wraplength=wraplength,
                               cex.comps=cex.comps,
                               ticklabels=ticklabels)
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Comparisons Graph, Resistant & Robust analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              comparisonsGraphStamp(mcadjust, alphapercent,
                                    desc=settings$type)
            }
            
            else if(aft && device=="single") {

              comparisonsgraph(aftcomprs, difftype, analysisname,
                               endptlabel, alpha, digits, titlestamp=FALSE,
                               wraplength=wraplength,
                               cex.comps=cex.comps,
                               explanation=FALSE, ticklabels=ticklabels)
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Comparisons Graph, AFT analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              comparisonsGraphStamp(mcadjust, alphapercent,
                                    desc=settings$type)
            }
            else if(uv && device=="single") {
              comparisonsgraph(uvcomprs, difftype, analysisname,
                               endptlabel, alpha, digits, titlestamp=FALSE,
                               explanation=FALSE,
                               wraplength=wraplength,
                               cex.comps=cex.comps,
                               ticklabels=ticklabels)
              if(stamps) graphStampCG(grid=FALSE)
              ## Text Annotations
              title(main=paste("Comparisons Graph, Unequal Variances analysis\n",
                      analysisname, sep=""), line=2, cex.main=1.1)
              comparisonsGraphStamp(mcadjust, alphapercent,
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





validErrorDf <- function(errordf, varcovmatrix,
                         n, mcadjust) {
  ## 
  if(errordf=="approx") {
    ## First check that estimates are uncorrelated, i.e
    ## the varcovmatrix is diagonal.
    if(sum(varcovmatrix[lower.tri(varcovmatrix)],
           varcovmatrix[upper.tri(varcovmatrix)])!=0) {
      stop(cgMessage("The varcovmatrix argument appears to have",
                     "nonzero off-diagonal elements. Currently,",
                     "the \"approx\" setting for the errordf",
                     "argument will only work for uncorrelated estimates."))
    }
    if(!isAllEqual(c(length(diag(varcovmatrix)),
                     length(n)))) {
      stop(cgMessage("The lengths of the n argument and",
                     "the number of diagonal elements in the varcovmatrix",
                     "argument appear to differ. They need to be the same."))
    }
    if(mcadjust) {
      stop(cgMessage("The mcadjust argument must be FALSE when the",
                     "errordf argument is \"approx\"."))
    }
  }
  else if(length(errordf)!=1 || !is.numeric(errordf) || errordf < 1) {
    stop(cgMessage("The errordf argument needs to be length 1",
                   "and set to be 1 or greater or the value \"approx\"."))
  }
  else return(TRUE)
}


validComparisonType <- function(type, errordf) {
  type <- try(match.arg(type, c("pairwisereflect","pairwise",
                                "allgroupstocontrol", "custom")))
  if(class(type)=="try-error") {
    stop(cgMessage("The type argument must be one of \"pairwisereflect\",",
                   "\"pairwise\", \"allgroupstocontrol\",",
                   "or \"custom\"."))
  }
  
  if(errordf=="approx" && !any(is.element(type,
       c("pairwisereflect","pairwise",
         "allgroupstocontrol")))) {
    stop(cgMessage("The type argument must be one of \"pairwisereflect\",",
                   "\"pairwise\", or \"allgroupstocontrol\"",
                   "when the errordf argument is \"approx\"."))
  }
  else if(!any(is.element(type, c("pairwisereflect","pairwise",
                                  "allgroupstocontrol", "custom")))) {
    stop(cgMessage("The type argument must be one of \"pairwisereflect\",",
                   "\"pairwise\", \"allgroupstocontrol\",",
                   "or \"custom\"."))
  }
  else return(TRUE)
  
}


validEstimates <- function(x) {
  if(!is.numeric(x) || length(x) < 2 ) {
    stop(cgMessage("The estimates argument needs to be numeric",
                   "and be of length two or more."
                   ))
  }
  if(is.null(thenames <- names(x))) {
    stop(cgMessage("The estimates arguments needs to have",
                   "a names attribute to use for",
                   "thes group labels."))
  }
  if(any(is.na(thenames)) || any(thenames=="")) {
    stop(cgMessage("One or more of the groups names in the",
                   "names attribute of estimates",
                   "appears to be missing."))
  }
  else return(x)
}



validGrpNames <- function(x) {
  if(!is.character(x) || length(x) < 2 ) {
    stop(cgMessage("The grpnames argument needs to be character",
                   "and be of length two or more.",
                   ))
  }
  if(any(is.na(x)) || any(x=="")) {
    stop(cgMessage("One or more of the groups names in the",
                   "names attribute of estimates",
                   "appears to be missing.",
                   ))
  }
  else return(x)
}



validN <- function(n, estimates) {
  if(!is.numeric(n) || any(n < 1)) {
    stop(cgMessage("The n argument vector needs to be numeric",
                   "and all elements must be one or greater."))
  }
  else if(length(n)!=length(estimates)) {
    stop(cgMessage("The n argument vector needs to be of the",
                   "same length as the estimates vector."))
  }
  else return(as.integer(round(n, 0)))
}



validCnames <- function(type, x, reqlength) {
  if(length(x)==1 && x=="derive") return(TRUE)
  if(type!="custom") {
    stop(cgMessage("The cnames argument must be left at its",
                   "default value of \"derive\" when type is",
                   "not \"custom\"."))
  }
  if(length(x) != reqlength) {
    stop(cgMessage("The cnames argument must be either the",
                   "value \"derive\", ",
                   "or must be of length", reqlength,
                   ", the number of rows of the",
                   "contrastmatrix argument."))
  }
  return(TRUE)
}

validAddPct <- function(addpct, endptscale) {
  if(addpct && endptscale=="log") {
    stop(cgMessage("The addpct argument cannot be TRUE,",
                   "and must be set to FALSE,",
                   "when the \"endptscale\" argument is",
                   "\"log\"."))

  }
  return(TRUE)
}



