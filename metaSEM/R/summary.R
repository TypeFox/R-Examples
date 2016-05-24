summary.wls <- function(object, df.adjustment=0, R=50, ...) {
    if (!is.element("wls", class(object)))
      stop("\"object\" must be an object of class \"wls\".")

    n <- object$n
    tT <- object$mx.fit@output$Minus2LogLikelihood
    my.mx <- summary(object$mx.fit)
    ## Adjust the df by the no. of constraints on the diagonals
    ## Adjust the df manually with df.adjustment
    if (object$diag.constraints) {
      dfT <- object$noObservedStat - my.mx$estimatedParameters + sum(object$Constraints) + df.adjustment
      no.constraints <- sum(object$Constraints)
    } else {
      dfT <- object$noObservedStat - my.mx$estimatedParameters + df.adjustment
      no.constraints <- 0
    }
    tB <- object$indepModelChisq
    dfB <- object$indepModelDf
    p <- pchisq(tT, df=dfT, lower.tail=FALSE)
    sampleS <- mxEval(sampleS, object$mx.fit)
    impliedS <- mxEval(impliedS, object$mx.fit)

    cor.analysis <- object$cor.analysis

    ## Hu and Bentler (1998) Psychological Methods
    ## Protect RMSEA divided by 0
    if (dfT==0) {
      RMSEA <- 0
    } else {
      RMSEA <- sqrt(max((tT-dfT)/(n-1),0)/dfT)
    }
    ## RMSEA 95% CI
    RMSEA.CI <- .rmseaCI(chi.squared=tT, df=dfT, N=n, G=1, lower=.05, upper=.95)
    
    TLI <- (tB/dfB - tT/dfT)/(tB/dfB-1)
    CFI <- 1 - max((tT-dfT),0)/max((tT-dfT),(tB-dfB),0)
	## FIXME: better to use -2LL+2r where r is the no. of free parameters (Mplus, p. 22; R ?AIC)
	## This will be more consistent with logLik(). However, it seems that df not npar is used in logLik().
    AIC <- tT-2*dfT
    BIC <- tT-log(n)*dfT
    if (cor.analysis) {
      # diag(impliedS)!=1
      SRMR <- sqrt(mean(vechs(sampleS-impliedS)^2))
    } else {
      # standardize it according to Hu & Bentler (1998)
      stand <- Diag(1/sqrt(Diag(sampleS)))
      SRMR <- sqrt(mean(vech(stand %*% (sampleS-impliedS) %*% stand)^2))
    }

    stat <- matrix(c(n, tT, dfT, p, no.constraints, df.adjustment, tB, dfB,
                     RMSEA, RMSEA.CI[1], RMSEA.CI[2], SRMR, TLI, CFI, AIC, BIC), ncol=1)

    dimnames(stat) <- list( c("Sample size", "Chi-square of target model", "DF of target model",
                        "p value of target model", "Number of constraints imposed on \"Smatrix\"",
                        "DF manually adjusted", "Chi-square of independence model",
                        "DF of independence model",  "RMSEA", "RMSEA lower 95% CI",
                        "RMSEA upper 95% CI", "SRMR", "TLI", "CFI", "AIC", "BIC"), "Value" )
   
    ## my.para <- my.mx$parameters       # Worked up to OpenMx1.0.6
    my.para <- my.mx$parameters[, 1:6]   # Fixed for OpenMx1.1
    
    ## ## Names in matrix, e.g., P[1,2], L[1,2], ...
    ## my.matrices <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    ## ## Names in labels
    ## my.labels <- my.para$name
    ## ## Remove "model.name". in my.labels
    ## my.labels <- gsub(paste(object$mx.model$name, ".", sep=""), "", my.labels)

    ## Check if CIs on parameter estimates are present
    ## Unsafe to check my.mx$CI
    if (object$intervals.type=="z") {
        
      ## Since no SE, parametric bootstrap is used
      if (object$diag.constraints)
        my.para$Std.Error <- sqrt(Diag(vcov.wls(object, R=R)))
        
      my.para$lbound <- my.para$Estimate - qnorm(.975)*my.para$Std.Error
      my.para$ubound <- my.para$Estimate + qnorm(.975)*my.para$Std.Error

      ## Order the parameters according to matrices, row and col
      my.para <- my.para[order(my.para$matrix, my.para$row, my.para$col), , drop=FALSE]
      
      ## Remove "model.name". in my.labels in case no labels are given 
      my.labels <- gsub(paste(object$mx.model$name, ".", sep=""), "", my.para$name)
      ## Report labels, not matrices
      dimnames(my.para)[[1]] <- my.labels     
      
      coefficients <- my.para[, -c(1:4)]      
    } else {
      ## # model.name: may vary in diff models
      ## model.name <- object$call[[match("model.name", names(object$call))]]
      ## # if not specified, the default is "Structure" as it can be either "Correlation Structure" or "Covariance Structure" 
      ## if (is.null(model.name)) {
      ##   model.name <- "Structure."
      ## } else {
      ##   model.name <- paste(model.name, ".", sep="")
      ## }          
      ## name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
      ##                  {strsplit(x, model.name, fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)

      ## Simply remove the part before "."
      ## name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
      ##                  {strsplit(x, ".", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)

      ## Convert a data frame with length of 0 in my.mx$CI and remove the last column "note"
      my.ci <- my.mx$CI
      if (length(my.ci)==0) my.ci <- NULL else my.ci <- my.ci[,1:3, drop=FALSE]        

      ## Labels, not matrices  
      name <- dimnames(my.ci)[[1]]
        
      my.ci <- data.frame(name, my.ci)
      my.para <- merge(my.para, my.ci, by=c("name"))
      
      ## Order the parameters according to matrices, row and col
      my.para <- my.para[order(my.para$matrix, my.para$row, my.para$col), , drop=FALSE]
      
      ## Remove "model.name". in my.labels in case no labels are given
      my.labels <- gsub(paste(object$mx.model$name, ".", sep=""), "", my.para$name)
      
      ## Report labels, not matrices
      dimnames(my.para)[[1]] <- my.labels
      
      # NA for LBCI
      my.para$Std.Error <- NA
      
      coefficients <- my.para[, -c(1:4, 8)] 
    }
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))
    
    ## Extract mx.algebras (and CI)
    if (is.null(object$mx.algebras)) {
      mx.algebras <- NULL
    } else {
      if (object$intervals.type=="z") {
        mx.algebras <- object$mx.fit@algebras[object$mx.algebras]
        mx.algebras <- sapply(mx.algebras, function(x) { if (!is.null(x)) x@result})
        ## remove lists with NA names
        mx.algebras <- unlist(mx.algebras[!is.na(names(mx.algebras))])
      } else {
          
        ## Convert a data frame with length of 0 in my.mx$CI and remove the last column "note"
        my.ci <- my.mx$CI
        if (length(my.ci)==0) my.ci <- NULL else my.ci <- my.ci[,1:3, drop=FALSE]
          
        ## Remove the first part of "model name"."algebra" in row names
        ## return NA is no "model name". in row names
        name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                        {strsplit(x, ".", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)

        # dimnames(my.ci)[[1]] <- name
        ## my.matches <- grep(paste(object$mx.algebras, 
        ##         "\\[", sep = "", collapse = "|"), row.names(my.ci), value = TRUE)
        # mx.algebras <- my.ci[!is.na(name), , drop = FALSE]
        
        ## Exclude the parameters and keep the functions
        mx.algebras <- my.ci[!is.na(name), , drop=FALSE]
        dimnames(mx.algebras)[[1]] <- name[!is.na(name)]
        
        dimnames(mx.algebras)[[2]] <- c("lbound", "Estimate", "ubound")
      }
    }

    out <- list(call=object$call, coefficients=coefficients, stat=stat, intervals.type=object$intervals.type,
                Mx.status1=object$mx.fit@output$status[[1]], mx.algebras=mx.algebras)
    class(out) <- "summary.wls"
    out
}

print.summary.wls <- function(x, ...) {
    if (!is.element("summary.wls", class(x)))
    stop("\"x\" must be an object of class \"summary.wls\".")
    ## call.text <- deparse(x$call)
    ## cat("Call:\n")
    ## for (i in 1:length(call.text)) {
    ##     cat(call.text[i], "\n")
    ## }	

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("95% confidence intervals: ")
    switch(x$intervals.type,
           z = cat("z statistic approximation"),
           LB = cat("Likelihood-based statistic") )

    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    if (!is.null(x$mx.algebras)) {
      if (x$intervals.type=="z") {
        cat("\nmxAlgebras objects:\n")
        print(x$mx.algebras)
      } else {
        cat("\nmxAlgebras objects (and their 95% likelihood-based CIs):\n")
        print(x$mx.algebras)
      }
    }
    
    cat("\nGoodness-of-fit indices:\n")
    printCoefmat(x$stat, ...)
    
    ## cat("\nR version:", x$R.version)
    ## cat("\nOpenMx version:", x$OpenMx.version)
    ## cat("\nmetaSEM version:", x$metaSEM.version)
    ## cat("\nDate of analysis:", x$date)
    cat("OpenMx status1:", x$Mx.status1, "(\"0\" or \"1\": The optimization is considered fine.\nOther values indicate problems.)\n")
    ## cat("OpenMx status1:", x$Mx.status1, "(\"0\" and \"1\": considered fine; other values indicate problems)\n")
    ## cat("\nSee http://openmx.psyc.virginia.edu/wiki/errors for the details.\n\n")
}

## NOTE:
## 1. The chi-square of the independence model in LISREL is different from other SEM programs.
## Thus, some GOF in stage 1 analysis are different from tssem1FEM()
## 2. When there are missing correlations or covariances in TSSEM program (LISREL),
## they are treated as observed correlations or covariances without any constraint.
## Thus, the DF of the independence model and some GOF in LISREL are incorrect.
summary.tssem1FEM <- function(object, ...) {
    if (!is.element("tssem1FEM", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM\".")
    
    cor.analysis <- object$cor.analysis
    
    pooledS <- coef.tssem1FEM(object)
    no.var <- ncol(pooledS)
    
    # Fixed if there are incomplete data in the first group
    ## no.var <- ncol(mxEval(S1, object$mx.fit))
    if (cor.analysis)
      ps <- no.var*(no.var-1)/2 else ps <- no.var*(no.var+1)/2

    ## A vector of sample sizes
    n <- object$n

    ## Check error in running independence model
    if (inherits(object$baseMinus2LL, "error")) {
      tT <- NA
      tB <- NA
      dfB <- NA
    } else {
      tT <- object$mx.fit@output$Minus2LogLikelihood - unname(object$baseMinus2LL["SaturatedLikelihood"])
      tB <- unname(object$baseMinus2LL["IndependenceLikelihood"]) - unname(object$baseMinus2LL["SaturatedLikelihood"])
      dfB <- unname(object$baseMinus2LL["independenceDoF"])
    }

    my.fit <- summary(object$mx.fit)
    ## tT <- object$modelMinus2LL - object$saturatedMinus2LL 
    dfT <- my.fit$degreesOfFreedom
    ## tB <- object$independentMinus2LL - object$saturatedMinus2LL
    ## dfB <- mx.fit$degreesOfFreedom + ps
    p <- pchisq(tT, df=dfT, lower.tail=FALSE)

    no.groups <- length(object$mx.fit@submodels)
    ## Steiger (1998, Eq. 24) A note on multiple sample extensions of the RMSEA fit indices. SEM, 5(4), 411-419.
    RMSEA <- sqrt(no.groups)*sqrt(max((tT-dfT)/(sum(n)-1),0)/dfT)
    ## RMSEA 95% CI
    RMSEA.CI <- .rmseaCI(chi.squared=tT, df=dfT, N=sum(n), G=no.groups, lower=.05, upper=.95)   
    
    ## Hu and Bentler (1998)
    TLI <- (tB/dfB - tT/dfT)/(tB/dfB-1)
    CFI <- 1 - max((tT-dfT),0)/max((tT-dfT),(tB-dfB),0)
    ## FIXME: definitions of AIC
    AIC <- tT-2*dfT
    BIC <- tT-log(sum(n))*dfT

    ## SRMR per study
    ## Index for missing variables: only check the diagonals only!!!
    #miss.index <- lapply(x$data, function(x) { is.na(diag(x)) })
    srmr <- function(sampleS, pooledS, cor.analysis) {
      index <- is.na(Diag(sampleS))
      # Selection matrix to select complete data
      Sel <- Diag(rep(1, ncol(pooledS)))
      Sel <- Sel[!index, ]
      sampleS <- sampleS[!index, !index]
      if (cor.analysis) {
        sqrt(mean( vechs( cov2cor(sampleS)- Sel %*% pooledS %*% t(Sel) )^2 ))
      } else {
        stand <- Diag(1/sqrt(Diag(sampleS)))
        sqrt(mean( vech( stand %*% (sampleS - Sel %*% pooledS %*% t(Sel) ) %*% stand )^2 ))
      }
    }

    ## weight to calculate SRMR
    n.weight <- (n-1)/sum(n-1)
    SRMR <- sapply(object$data, srmr, pooledS=pooledS, cor.analysis=cor.analysis)
    SRMR <- sum(n.weight*SRMR)
    ## SRMR <- sqrt(mean(unlist(sapply(object$data, srmr, pooledS=object$pooledS, cor.analysis=cor.analysis))))
    #SRMR <- sqrt(mean(mapply(srmr, x$data, miss.index,
    #                         MoreArgs=list(pooledS=x$pooledS, cor.analysis=cor.analysis))))
    
    stat <- matrix(c(sum(n), tT, dfT, p, tB, dfB, RMSEA, RMSEA.CI[1], RMSEA.CI[2],
                     SRMR, TLI, CFI, AIC, BIC), ncol=1)
    rownames(stat) <- c("Sample size", "Chi-square of target model", "DF of target model",
                        "p value of target model", "Chi-square of independence model",
                        "DF of independence model", "RMSEA", "RMSEA lower 95% CI",
                        "RMSEA upper 95% CI","SRMR", "TLI", "CFI", "AIC", "BIC")
    colnames(stat) <- "Value"

    # calculate coefficients    
    my.para <- my.fit$parameters
    ## my.para <- my.para[my.para$matrix=="S1", ]
    my.para <- my.para[my.para$matrix=="S", , drop=FALSE]
    #Sel <- grep("^S", my.para$matrix, value=TRUE)
    #my.para <- subset(my.para, my.para$matrix==Sel)
    my.para <- my.para[order(my.para$row, my.para$col), , drop=FALSE]
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    dimnames(my.para)[[1]] <- my.para$name
    coefficients <- my.para[, c(5,6)]
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))
    
    out <- list(call=object$call, coefficients=coefficients, stat=stat, Mx.status1=object$mx.fit@output$status[[1]])
    class(out) <- "summary.tssem1FEM"
    out
}
   
print.summary.tssem1FEM <- function(x, ...) {
    if (!is.element("summary.tssem1FEM", class(x)))
    stop("\"x\" must be an object of class \"summary.tssem1FEM\".")
    ## call.text <- deparse(x$call)
       
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("Coefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    cat("\nGoodness-of-fit indices:\n")
    printCoefmat(x$stat, ...)

    ## cat("\nR version:", x$R.version)
    ## cat("\nOpenMx version:", x$OpenMx.version)
    ## cat("\nmetaSEM version:", x$metaSEM.version)
    ## cat("\nDate of analysis:", x$date)
    cat("OpenMx status1:", x$Mx.status1, "(\"0\" or \"1\": The optimization is considered fine.\nOther values may indicate problems.)\n")
    ## cat("OpenMx status:", x$Mx.status1, "(\"0\" and \"1\": considered fine; other values indicate problems)\n")
    ## cat("\nSee http://openmx.psyc.virginia.edu/wiki/errors for the details.\n\n")    
}

print.tssem1FEM <- function(x, ...) {
    if (!is.element("tssem1FEM", class(x)))
      stop("\"x\" must be an object of class \"tssem1FEM\".")
    ## call.text <- deparse(x$call)

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("Structure:\n")
    print(summary.default(x), ...)
}

print.tssem1FEM.cluster <- function(x, ...) {
    if (!is.element("tssem1FEM.cluster", class(x)))
      stop("\"x\" must be an object of class \"tssem1FEM.cluster\".")

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("Structure:\n")
    print(summary.default(x), ...)
}

print.tssem1REM <- function(x, ...) {
    if (!is.element("tssem1REM", class(x)))
      stop("\"x\" must be an object of class \"tssem1REM\".")
    ## call.text <- deparse(x$call)

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("Structure:\n")
    print(summary.default(x), ...)
}

summary.tssem1REM <- function(object, ...) {
  if (!is.element("tssem1REM", class(object)))
    stop("\"object\" must be an object of class \"tssem1REM\".")
  summary.meta(object)
}
   
print.wls <- function(x, ...) {
    if (!is.element("wls", class(x)))
      stop("\"x\" must be an object of class \"wls\".")

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("Structure:\n")	
    print(summary.default(x), ...)
}

print.meta <- function(x, ...) {
    if (!is.element("meta", class(x)))
      stop("\"x\" must be an object of class \"meta\".")
    ## cat("Call:\n")
    ## cat(deparse(x$call))
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")     
    cat("Structure:\n")
    print(summary.default(x), ...)
}

summary.meta <- function(object, homoStat=TRUE, ...) {
    if (!is.element("meta", class(object)))
      stop("\"object\" must be an object of class \"meta\".")

    # calculate coefficients    
    my.mx <- summary(object$mx.fit)
    ## Exclude lbound ubound etc
    my.para <- my.mx$parameters[, 1:6, drop=FALSE]   
    ## OpenMx1.1: y1, y2 and x1 appear in col
    ## my.para$col <- sub("[a-z]", "", my.para$col)  # Fixed for OpenMx1.1
    ## For example, P[1,2], L[1,2], ...
    ## my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    
    # Determine if CIs on parameter estimates are present
    if (object$intervals.type=="z") {
      my.para$lbound <- my.para$Estimate - qnorm(.975)*my.para$Std.Error
      my.para$ubound <- my.para$Estimate + qnorm(.975)*my.para$Std.Error
      # remove rows with missing labels
      ## my.para <- my.para[!is.na(my.para$label), ,drop=FALSE]
      ## my.para <- my.para[order(my.para$label), , drop=FALSE]
      coefficients <- my.para[, -c(1:4), drop=FALSE]
      dimnames(coefficients)[[1]] <- my.para$name
    } else {
      ## # model.name: may vary in diff models
      ## ## FIXME: it may break down when model.name is a variable.
      ## model.name <- object$call[[match("model.name", names(object$call))]]
      ## # if not specified, the default in meta() is "Meta analysis with ML"
      ## if (is.null(model.name)) {
      ##   model.name <- "Meta analysis with ML."
      ## } else {
      ##   model.name <- paste(model.name, ".", sep="")
      ## }
      ## name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
      ##                  {strsplit(x, model.name, fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
        
      ## name <- sapply(unlist(dimnames(my.ci)[1]), function(x) {strsplit(x, ".", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)

      ## Convert a data frame with length of 0 in my.mx$CI and remove the last column "note"
      my.ci <- my.mx$CI
      if (length(my.ci)==0) my.ci <- NULL else my.ci <- my.ci[,1:3, drop=FALSE]        
        
      ## Select the elements matched my.para (excluded I2)  
      my.ci <- my.ci[row.names(my.ci) %in% my.para$name, ]
      
      my.ci <- data.frame(name=row.names(my.ci), my.ci)
      my.para <- merge(my.para, my.ci, by=c("name"))      
      ## # remove rows with missing labels
      ## my.para <- my.para[!is.na(my.para$label), , drop=FALSE]
      ## my.para <- my.para[order(my.para$label), , drop=FALSE]
      coefficients <- my.para[, -c(1:4,8)]
      dimnames(coefficients)[[1]] <- my.para$name
      # NA for LBCI
      coefficients$Std.Error <- NA
      ## coefficients$"z value" <- NA
      ## coefficients$"Pr(>|z|)" <- NA         
    }
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))

    intervals.type <- object$call[[match("intervals.type", names(object$call))]]
    # default
    if (is.null(intervals.type))
      intervals.type <- "z"

    ## Homogeneity statistic
    no.y <- object$no.y
    no.v <- no.y*(no.y+1)/2
	# Remove studies that have missing x. Make sure that studies are the same in calculating Q.stat and meta()
    if (homoStat) {
      Q.stat <- homoStat(y=object$data[!object$miss.x, 1:no.y], v=object$data[!object$miss.x, (no.y+1):(no.y+no.v)])
    } else {
      Q.stat <- list(Q=NA, Q.df=NA, pval=NA)
    }

    ## Calculate I2
    if (object$no.x==0) {
      ## I2.values <- .I2(object, my.mx, my.ci, model.name)
      I2.values <- .I2(object, my.mx)
      R2.values <- NA
    } else {## no.x != 0
      I2.values <- NA
      if (object$R2) {
        R2.values <- .R2(object)
        ## Ensure R2 is within 0 and 1
        R2.values <- ifelse(R2.values < 0, 0, R2.values)
        R2.values <- ifelse(R2.values > 1, 1, R2.values)
      } else {
        R2.values <- NA
      }
    }
          
    out <- list(call=object$call, Q.stat=Q.stat, intervals.type=intervals.type,
                I2=object$I2, I2.values=I2.values, R2=object$R2, R2.values=R2.values, no.studies=my.mx$numObs,
                obsStat=my.mx$observedStatistics, estPara=my.mx$estimatedParameters,
                df=my.mx$degreesOfFreedom, Minus2LL=my.mx$Minus2LogLikelihood,
                coefficients=coefficients, Mx.status1=object$mx.fit@output$status[[1]])
    class(out) <- "summary.meta"
    out
}

## Magee, L. (1990). R2 Measures Based on Wald and Likelihood Ratio Joint Significance Tests. The American Statistician, 44(3), 250-253. doi:10.2307/2685352
## Nagelkerke, N. J. D. (1991). A note on a general definition of the coefficient of determination. Biometrika, 78(3), 691-692. doi:10.2307/2337038
## Minus2LLbase <- summary(object$mx0.fit$mx.fit)$Minus2LogLikelihood 
## Minus2LLmodel <- my.mx$Minus2LogLikelihood           
## Minus2LLbase, Minus2LLmodel, (1 - exp((Minus2LLmodel-Minus2LLbase)/no.studies))
## "-2LL (no predictor)", "-2LL (with predictors)", "R2 (pseudo)"

print.summary.meta <- function(x, ...) {
    if (!is.element("summary.meta", class(x)))
    stop("\"x\" must be an object of class \"summary.meta\".")

    ## cat("Call:\n")
    ## cat(deparse(x$call))
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    
    cat("95% confidence intervals: ")
    switch(x$intervals.type,
           z = cat("z statistic approximation"),
           LB = cat("Likelihood-based statistic") )

    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    cat("\nQ statistic on the homogeneity of effect sizes:", x$Q.stat[["Q"]])
    cat("\nDegrees of freedom of the Q statistic:", x$Q.stat[["Q.df"]])
    cat("\nP value of the Q statistic:", x$Q.stat[["pval"]])

    ## Print heterogeneity indices if no x in the call
    if ( is.null(x$call[[match("x", names(x$call))]]) ) {
      switch(x$intervals.type,
             z =  { cat("\n\nHeterogeneity indices (based on the estimated Tau2):\n")
                    printCoefmat(x$I2.values[,2, drop=FALSE], ...) },
             LB = { cat("\n\nHeterogeneity indices (I2) and their 95% likelihood-based CIs:\n")
                    printCoefmat(x$I2.values, ...) } )
    } else {
      ## There are predictors
      if (x$R2) {
        cat("\n\nExplained variances (R2):\n")
        printCoefmat(x$R2.values, ...) 
      }
    }    
    cat("\nNumber of studies (or clusters):", x$no.studies)
    cat("\nNumber of observed statistics:", x$obsStat)
    cat("\nNumber of estimated parameters:", x$estPara)
    cat("\nDegrees of freedom:", x$df)
    cat("\n-2 log likelihood:", x$Minus2LL, "\n")        
    cat("OpenMx status1:", x$Mx.status1, "(\"0\" or \"1\": The optimization is considered fine.\nOther values may indicate problems.)\n")
    ## cat("OpenMx status1:", x$Mx.status1, "(\"0\" and \"1\": considered fine; other values indicate problems)\n")
    ## cat("\nSee http://openmx.psyc.virginia.edu/wiki/errors for the details.\n\n")      
}


## summary.reml <- function(object, ...) {
##     if (!is.element("reml", class(object)))
##     stop("\"object\" must be an object of class \"reml\".")

##     # calculate coefficients    
##     my.mx <- summary(object$mx.fit)
##     ## my.para <- my.mx$parameters       # Worked up to OpenMx1.0.6
##     my.para <- my.mx$parameters[, 1:6]   # Fixed for OpenMx1.1 
##     # For example, P[1,2], L[1,2], ...
##     my.para$label <- my.para$name
##     my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    
##     my.ci <- my.mx$CI
##     # Determine if CIs on parameter estimates are present
##     if (is.null(dimnames(my.ci))) {
##       my.para$lbound <- my.para$Estimate - qnorm(.975)*my.para$Std.Error
##       my.para$ubound <- my.para$Estimate + qnorm(.975)*my.para$Std.Error
##       my.para <- my.para[order(my.para$matrix, my.para$row, my.para$col), ]
##       # remove rows with missing labels
##       my.para <- my.para[!is.na(my.para$label), ]
##       coefficients <- my.para[, -c(1:4,7)]
##       dimnames(coefficients)[[1]] <- my.para$label 
##     } else {
##       # model.name: may vary in diff models
##       model.name <- object$call[[match("model.name", names(object$call))]]
##       # if not specified, the default in meta() is "Variance component with REML"
##       if (is.null(model.name)) {
##         model.name <- "Variance component with REML."
##       } else {
##         model.name <- paste(model.name, ".", sep="")
##       }      
##       name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
##                        {strsplit(x, model.name, fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
##       # remove duplicate elements in my.ci from my.para$name
##       name.sel <- name %in% my.para$name
##       my.ci <- data.frame(name=my.para$name, my.ci[name.sel, ,drop=FALSE])
##       my.para <- merge(my.para, my.ci, by=c("name"))
##       my.para <- my.para[order(my.para$matrix, my.para$row, my.para$col), ]
##       # remove rows with missing labels
##       my.para <- my.para[!is.na(my.para$label), ]
##       coefficients <- my.para[, -c(1:4,7,9)]
##       dimnames(coefficients)[[1]] <- my.para$label 
##     }
##     coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
##     coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))

##     intervals.type <- object$call[[match("intervals.type", names(object$call))]]
##     # default
##     if (is.null(intervals.type))
##       intervals.type <- "z"

##     out <- list(call=object$call, intervals.type=intervals.type, no.studies=my.mx$numObs,
##                 obsStat=my.mx$observedStatistics, estPara=my.mx$estimatedParameters,
##                 df=my.mx$degreesOfFreedom, Minus2LL=my.mx$Minus2LogLikelihood,
##                 coefficients=coefficients, Mx.status1=object$mx.fit@output$status[[1]])
##     class(out) <- "summary.reml"
##     out
## }

summary.reml <- function(object, ...) {
    if (!is.element("reml", class(object)))
    stop("\"object\" must be an object of class \"reml\".")

    # calculate coefficients    
    my.mx <- summary(object$mx.fit)
    ## my.para <- my.mx$parameters       # Worked up to OpenMx1.0.6
    my.para <- my.mx$parameters[, 1:6]   # Fixed for OpenMx1.1 
    ## # For example, P[1,2], L[1,2], ...
    ## my.para$label <- my.para$name   
    
    # Determine if CIs on parameter estimates are present
    if (object$intervals.type=="z") {
      my.para$lbound <- my.para$Estimate - qnorm(.975)*my.para$Std.Error
      my.para$ubound <- my.para$Estimate + qnorm(.975)*my.para$Std.Error      
      # remove rows with missing labels
      ## my.para <- my.para[!is.na(my.para$label), ]
      ## my.para <- my.para[order(my.para$matrix, my.para$row, my.para$col), ]
      coefficients <- my.para[, -c(1:4)]
      dimnames(coefficients)[[1]] <- my.para$name
    } else {

      ## if ( "reml3" %in% class(object) ) {
      ##   name.sel <- my.para$name
      ## } else {
      ##   my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
      ##   ## # model.name: may vary in diff models
      ##   ## model.name <- object$call[[match("model.name", names(object$call))]]
      ##   ## # if not specified, the default in meta() is "Variance component with REML"
      ##   ## if (is.null(model.name)) {
      ##   ##   model.name <- "Variance component with REML."
      ##   ## } else { # reml2
      ##   ##   model.name <- paste(model.name, ".", sep="")
      ##   ## }      
      ##   ## name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
      ##   ##                  {strsplit(x, model.name, fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
      ##   name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
      ##                        {strsplit(x, ".", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
        
      ##   # remove duplicate elements in my.ci from my.para$name
      ##   name.sel <- name %in% my.para$name
      ## } # reml2

      ## Convert a data frame with length of 0 in my.mx$CI and remove the last column "note"
      my.ci <- my.mx$CI
      if (length(my.ci)==0) my.ci <- NULL else my.ci <- my.ci[,1:3, drop=FALSE]        

      my.ci <- data.frame(name=row.names(my.ci), my.ci)
      my.para <- merge(my.para, my.ci, by=c("name"))      
      ## # remove rows with missing labels
      ## my.para <- my.para[!is.na(my.para$label), , drop=FALSE]
      ## my.para <- my.para[order(my.para$label), , drop=FALSE]
      coefficients <- my.para[, -c(1:4,8)]
      dimnames(coefficients)[[1]] <- my.para$name
      # NA for LBCI
      coefficients$Std.Error <- NA
    } 
      
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))

    ## intervals.type <- object$call[[match("intervals.type", names(object$call))]]
    ## # default
    ## if (is.null(intervals.type))
    ##   intervals.type <- "z"

    obsStat=object$numStats
    estPara=my.mx$estimatedParameters
    df=obsStat-estPara
    out <- list(call=object$call, intervals.type=object$intervals.type, no.studies=object$numObs,
                obsStat=obsStat, estPara=estPara,
                df=df, Minus2LL=my.mx$Minus2LogLikelihood,
                coefficients=coefficients, Mx.status1=object$mx.fit@output$status[[1]])
    class(out) <- "summary.reml"
    out
}

print.summary.reml <- function(x, ...) {
    if (!is.element("summary.reml", class(x)))
    stop("\"x\" must be an object of class \"summary.reml\".")
    
    ## cat("Call:\n")
    ## cat(deparse(x$call))
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("95% confidence intervals: ")
    switch(x$intervals.type,
           z = cat("z statistic approximation"),
           LB = cat("Likelihood-based statistic") )

    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    cat("\nNumber of studies (or clusters):", x$no.studies)
    cat("\nNumber of observed statistics:", x$obsStat)
    cat("\nNumber of estimated parameters:", x$estPara)
    cat("\nDegrees of freedom:", x$df)
    cat("\n-2 log likelihood:", x$Minus2LL, "\n")        
    cat("OpenMx status1:", x$Mx.status1, "(\"0\" or \"1\": The optimization is considered fine.\nOther values may indicate problems.)\n")
    ## cat("OpenMx status:", x$Mx.status1, "(\"0\" and \"1\": considered fine; other values indicate problems)\n")
    ## cat("\nSee http://openmx.psyc.virginia.edu/wiki/errors for the details.\n\n")      
}

print.reml <- function(x, ...) {
    if (!is.element("reml", class(x)))
      stop("\"x\" must be an object of class \"reml\".")
    ## cat("Call:\n")
    ## cat(deparse(x$call))

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")     
    cat("Structure:\n")
    print(summary.default(x), ...)
}

vcov.meta <- function(object, select=c("all", "fixed", "random"), ...) {
    if (!is.element("meta", class(object)))
    stop("\"object\" must be an object of class \"meta\".")

    # labels of the parameters    
    ## my.name <- summary(object$mx.fit)$parameters$name
    my.name <- names( omxGetParameters(object$mx.fit) )
    my.name <- my.name[!is.na(my.name)]
    
    select <- match.arg(select)
    switch( select,
         ## all = my.name <- my.name,
         fixed =  my.name <- my.name[ grep("Intercept|Slope", my.name) ],
         random = my.name <- my.name[ grep("Tau2", my.name) ]
         )
    
    # Fixed a bug that all elements have to be inverted before selecting some of them
    acov <- tryCatch( 2*solve(object$mx.fit@output$calculatedHessian)[my.name, my.name, drop=FALSE], error = function(e) e)
    # Issue a warning instead of error message
    if (inherits(acov, "error")) {
      cat("Error in solving the Hessian matrix.\n")
      warning(print(acov))
    } else {
      return(acov)
    }
}

vcov.tssem1FEM <- function(object, ...) {
  if (!is.element("tssem1FEM", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM\".")

    no.var <- ncol(coef.tssem1FEM(object))

    ## matrix of labels; only use the lower triangle
    ps.labels <- outer(1:no.var, 1:no.var, function(x, y) paste0("s", x, y))
    
    if (object$cor.analysis) {
        #Hessian_S <- 0.5*mx.fit@output$calculatedHessian[vechs(ps.labels), vechs(ps.labels)]
        acovS <- tryCatch( 2*solve(object$mx.fit@output$calculatedHessian)[vechs(ps.labels),
                           vechs(ps.labels)], error = function(e) e)
    } else {
        #Hessian_S <- 0.5*mx.fit@output$calculatedHessian[vech(ps.labels), vech(ps.labels)]
        acovS <-  tryCatch( 2*solve(object$mx.fit@output$calculatedHessian)[vech(ps.labels),
                            vech(ps.labels)], error = function(e) e)
    }
    # Issue a warning instead of error message
    if (inherits(acovS, "error")) {
      cat("Error in solving the Hessian matrix.\n")
      warning(print(acovS))
    } else {
      # Fixed a bug in a few lines later in dimnames(acovS) when acovS is a scalar
      acovS <- as.matrix(acovS)
    }

    acovS.dim <- outer(object$original.names, object$original.names, paste, sep="_")

    # create matrix of labels for ps
    if (object$cor.analysis) {
        psMatnames <- vechs(acovS.dim)
    } else {
        psMatnames <- vech(acovS.dim)
    }

    dimnames(acovS) <- list(psMatnames, psMatnames)

    acovS
}

vcov.tssem1FEM.cluster <- function(object, ...) {
  if (!is.element("tssem1FEM.cluster", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM.cluster\".")
  lapply(object, vcov.tssem1FEM)
}  

vcov.tssem1REM <- function(object, select=c("all", "fixed", "random"), ...) {
  if (!is.element("tssem1REM", class(object)))
    stop("\"object\" must be an object of class \"tssem1REM\".")
  vcov.meta(object, select, ...)
}

vcov.wls <- function(object, R=50, ...) {
  if (!is.element("wls", class(object)))
    stop("\"object\" must be an object of class \"wls\".")

  if (object$diag.constraints) {
    ## Parametric bootstrap for vcov when there are nonlinear constraints
    paraboot <- function(x) {
      if (x$cor.analysis) {
        sampleS <- mvrnorm(n=1, mu=vechs(x$Cov), Sigma=x$asyCov)
      } else {
        sampleS <- mvrnorm(n=1, mu=vech(x$Cov), Sigma=x$asyCov)
      }
      sampleS <- vec2symMat(sampleS, diag=!x$cor.analysis)
      mx.model <- mxModel(x$mx.fit, sampleS <- as.mxMatrix(sampleS, name="sampleS"))
      mx.model <- mxOption(mx.model, "Calculate Hessian", "No") 
      mx.model <- mxOption(mx.model, "Standard Errors"  , "No")  
      mx.fit <- tryCatch(mxRun(mx.model, intervals=FALSE, silent=TRUE, suppressWarnings=TRUE),
                         error = function(e) e )
      if (inherits(mx.fit, "error")) {
        return(NA)
      } else {
        omxGetParameters(mx.fit)
      }
    }
    #var(t(omxSapply(1:b, paraboot, simplify = TRUE)), na.rm=TRUE)
    acovS <- var(t(replicate(R, paraboot(object))), na.rm=TRUE)
    
    warning("Parametric bootstrap with ",R," replications was used to approximate the sampling covariance matrix of the parameter estimates. A better approach is to use likelihood-based confidence interval by including the intervals.type=\"LB\" argument in the analysis.\n")
        
  } else {
    ## Select the free parameters for inversion
    acovS <- tryCatch( 2*solve(object$mx.fit@output$calculatedHessian), error = function(e) e ) 

    # Issue a warning instead of error message
    if (inherits(acovS, "error")) {
      cat("Error in solving the Hessian matrix.\n")
      stop(print(acovS))
    } 
  }

  my.para <- summary(object$mx.fit)$parameters[, 1:4]
  my.labels <- my.para$name
  my.order <- with(my.para, order(matrix, row, col))

  acovS <- acovS[my.labels[my.order], my.labels[my.order]]
  my.labels <- my.labels <- gsub(paste(object$mx.model$name, ".", sep=""), "", row.names(acovS))
  dimnames(acovS) <- list(my.labels, my.labels)
  acovS
}

vcov.wls.cluster <- function(object, R=50, ...) {
    if (!is.element("wls.cluster", class(object)))
    stop("\"object\" must be an object of class \"wls.cluster\".")
    lapply(object, vcov.wls, R=R, ...)
}
  
vcov.reml <- function(object, ...) {
    if (!is.element("reml", class(object)))
    stop("\"object\" must be an object of class \"reml\".")

    ## # labels of the parameters    
    ## my.name <- summary(object$mx.fit)$parameters$name
    my.name <- names(omxGetParameters(object$mx.fit))
    my.name <- my.name[!is.na(my.name)]
    # Fixed a bug that all elements have to be inverted before selecting some of them
    acov <- tryCatch( 2*solve(object$reml@output$calculatedHessian)[my.name, my.name, drop=FALSE], error = function(e) e) 
    
    if (inherits(acov, "error")) {
      cat("Error in solving the Hessian matrix.\n")
      warning(print(acov))
    } else {
      return(acov)
    }
}
  

coef.meta <- function(object, select=c("all", "fixed", "random"), ...) {
  if (!is.element("meta", class(object)))
    stop("\"object\" must be an object of class \"meta\".")

  ## # labels of the parameters    
  ## my.name <- summary(object$mx.fit)$parameters$name
  ## my.name <- my.name[!is.na(my.name)]

  my.para <- omxGetParameters(object$mx.fit)
  select <- match.arg(select)
  switch( select,
         fixed =  my.para <- my.para[ grep("Intercept|Slope", names(my.para)) ],
         random = my.para <- my.para[ grep("Tau2", names(my.para)) ]
         )
  my.para
}

coef.tssem1FEM <- function(object, ...) {  
    if (!is.element("tssem1FEM", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM\".")

    pooledS <- eval(parse(text = "mxEval(S, object$mx.fit)"))
    dimnames(pooledS) <- list(object$original.names, object$original.names)  
    pooledS
}

coef.tssem1REM <- function(object, select=c("all", "fixed", "random"), ...) {
  if (!is.element("tssem1REM", class(object)))
    stop("\"object\" must be an object of class \"tssem1REM\".")

  coef.meta(object, select, ...)
}

coef.tssem1FEM.cluster <- function(object, ...) {
    if (!is.element("tssem1FEM.cluster", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM.cluster\".")
    lapply(object, coef.tssem1FEM)
}
  
coef.wls <- function(object, ...) {
    if (!is.element("wls", class(object)))
    stop("\"object\" must be an object of class \"wls\".")
    ## object$mx.fit@output$estimate
    my.mx <- summary(object$mx.fit)
    my.coef <- my.mx$parameters$Estimate
    ## # For example, P[1,2], L[1,2], ...
    ## names(my.coef) <- with(my.mx$parameters[, 2:4], paste(matrix,"[",row,",",col,"]",sep=""))
    my.labels <- my.mx$parameters$name
    my.labels <- gsub(paste(object$mx.model$name, ".", sep=""), "", my.labels)
    names(my.coef) <- my.labels    
    my.coef
}

coef.wls.cluster <- function(object, ...) {
    if (!is.element("wls.cluster", class(object)))
    stop("\"object\" must be an object of class \"wls.cluster\".")
    lapply(object, coef.wls)
}

coef.reml <- function(object, ...) {
    if (!is.element("reml", class(object)))
    stop("\"object\" must be an object of class \"reml\".")
    # labels of the parameters    
    ## my.name <- summary(object$mx.fit)$parameters$name
    my.name <- names(omxGetParameters(object$mx.fit))
    my.name <- my.name[!is.na(my.name)]
    object$mx.fit@output$estimate[my.name]
}

anova.meta <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$mx.fit)
  comparison <- lapply(list(...), function(x) x$mx.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

anova.wls <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$mx.fit)
  comparison <- lapply(list(...), function(x) x$mx.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

anova.reml <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$mx.fit)
  comparison <- lapply(list(...), function(x) x$mx.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

summary.tssem1FEM.cluster <- function(object, ...) {
    if (!is.element("tssem1FEM.cluster", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM.cluster\".")
    lapply(object, summary.tssem1FEM, ...)
}

summary.wls.cluster <- function(object, df.adjustment=0, R=50, ...) {
    if (!is.element("wls.cluster", class(object)))
    stop("\"object\" must be an object of class \"wls.cluster\".")
    lapply(object, summary.wls, df.adjustment=0, R=R, ...)
}

summary.meta3X <- function(object, allX=FALSE, ...) {
    if (!is.element("meta3X", class(object)))
      stop("\"object\" must be an object of class \"meta3X\".")

    # calculate coefficients    
    my.mx <- summary(object$mx.fit)
    ## Exclude lbound ubound etc
    my.para <- my.mx$parameters[, 1:6, drop=FALSE]   
    ## OpenMx1.1: y1, y2 and x1 appear in col
    ## my.para$col <- sub("[a-z]", "", my.para$col)  # Fixed for OpenMx1.1
    # For example, P[1,2], L[1,2], ...
    my.para$label <- my.para$name
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))

    if (!allX) {
      index.label <- my.para$label[grep("Intercept|Slope|Tau2", my.para$label)]
      index <- my.para$label %in% index.label
      my.para <- my.para[index, ]
    }

    ## Convert a data frame with length of 0 in my.mx$CI and remove the last column "note"
    my.ci <- my.mx$CI
    if (length(my.ci)==0) my.ci <- NULL else my.ci <- my.ci[,1:3, drop=FALSE]    
    
    # Determine if CIs on parameter estimates are present
    if (is.null(dimnames(my.ci))) {
      my.para$lbound <- my.para$Estimate - qnorm(.975)*my.para$Std.Error
      my.para$ubound <- my.para$Estimate + qnorm(.975)*my.para$Std.Error
      my.para <- my.para[order(my.para$label), , drop=FALSE]
      # remove rows with missing labels
      #my.para <- my.para[!is.na(my.para$label), ,drop=FALSE]
      coefficients <- my.para[, c("Estimate","Std.Error","lbound","ubound"), drop=FALSE]
    } else {
      name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                            {strsplit(x, ".", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
      my.ci <- data.frame(name, my.ci)
      my.para <- merge(x=my.para, y=my.ci, by=c("name"), all.x=TRUE)
      my.para <- my.para[order(my.para$label), , drop=FALSE]
      # remove rows with missing labels
      #my.para <- my.para[!is.na(my.para$label), , drop=FALSE]
      coefficients <- cbind(Estimate=my.para[, "Estimate"], Std.Error=NA, my.para[, c("lbound","ubound")])
      # NA for LBCI
      ## coefficients$Std.Error <- NA                             
    }
    dimnames(coefficients)[[1]] <- my.para$label
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))

    intervals.type <- object$call[[match("intervals.type", names(object$call))]]
    # default
    if (is.null(intervals.type))
      intervals.type <- "z"

    if (object$R2) {
      R2.values <- .R2(object)
    } else {
      R2.values <- NA
    }
              
    out <- list(call=object$call, intervals.type=intervals.type,
                R2=object$R2, R2.values=R2.values, no.studies=my.mx$numObs,
                obsStat=my.mx$observedStatistics, estPara=my.mx$estimatedParameters,
                df=my.mx$degreesOfFreedom, Minus2LL=my.mx$Minus2LogLikelihood,
                coefficients=coefficients, Mx.status1=object$mx.fit@output$status[[1]])
    class(out) <- "summary.meta3X"
    out
}

print.summary.meta3X <- function(x, ...) {
    if (!is.element("summary.meta3X", class(x)))
    stop("\"x\" must be an object of class \"summary.meta3X\".")

    ## cat("Call:\n")
    ## cat(deparse(x$call))
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    
    cat("95% confidence intervals: ")
    switch(x$intervals.type,
           z = cat("z statistic approximation"),
           LB = cat("Likelihood-based statistic") )

    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    if (x$R2) {
      cat("\nExplained variances (R2):\n")
      printCoefmat(x$R2.values, ...) 
    }
     
    cat("\nNumber of studies (or clusters):", x$no.studies)
    cat("\nNumber of observed statistics:", x$obsStat)
    cat("\nNumber of estimated parameters:", x$estPara)
    cat("\nDegrees of freedom:", x$df)
    cat("\n-2 log likelihood:", x$Minus2LL, "\n")        
    cat("OpenMx status1:", x$Mx.status1, "(\"0\" or \"1\": The optimization is considered fine.\nOther values may indicate problems.)\n")    
    ## cat("OpenMx status1:", x$Mx.status1, "(\"0\" and \"1\": considered fine; other values indicate problems)\n")
    ## cat("\nSee http://openmx.psyc.virginia.edu/wiki/errors for the details.\n\n")      
}

print.meta3X <- function(x, ...) {
    if (!is.element("meta3X", class(x)))
      stop("\"x\" must be an object of class \"meta3X\".")
    ## cat("Call:\n")
    ## cat(deparse(x$call))
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")     
    cat("Structure:\n")
    print(summary.default(x), ...)
}

vcov.meta3X <- function(object, select=c("all", "fixed", "random", "allX"), ...) {
    if (!is.element("meta3X", class(object)))
    stop("\"object\" must be an object of class \"meta3X\".")

    # labels of the parameters    
    ## my.name <- summary(object$mx.fit)$parameters$name
    my.name <- names( omxGetParameters(object$mx.fit) )
    my.name <- my.name[!is.na(my.name)]
    
    select <- match.arg(select)
    switch( select,
         all = my.name <- my.name[ grep("Intercept|Slope|Tau2", my.name) ],
         fixed =  my.name <- my.name[ grep("Intercept|Slope", my.name) ],
         random = my.name <- my.name[ grep("Tau2", my.name) ]
         )
    
    # Fixed a bug that all elements have to be inverted before selecting some of them
    acov <- tryCatch( 2*solve(object$mx.fit@output$calculatedHessian)[my.name, my.name, drop=FALSE], error = function(e) e)
    # Issue a warning instead of error message
    if (inherits(acov, "error")) {
      cat("Error in solving the Hessian matrix.\n")
      warning(print(acov))
    } else {
      return(acov)
    }
}

coef.meta3X <- function(object, select=c("all", "fixed", "random", "allX"), ...) {
  if (!is.element("meta3X", class(object)))
    stop("\"object\" must be an object of class \"meta3X\".")

  my.para <- omxGetParameters(object$mx.fit)
  select <- match.arg(select)
  switch( select,
         all =  my.para <- my.para[ grep("Intercept|Slope|Tau2", names(my.para)) ],
         fixed =  my.para <- my.para[ grep("Intercept|Slope", names(my.para)) ],
         random = my.para <- my.para[ grep("Tau2", names(my.para)) ]
         )
  my.para
}

anova.meta3X <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$mx.fit)
  comparison <- lapply(list(...), function(x) x$mx.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

coef.MxRAMModel <- function(object, ...) {
  if (!is.element("MxRAMModel", class(object)))
    stop("\"object\" must be an object of class \"MxRAMModel\".")
  omxGetParameters(object, ...)
}

vcov.MxRAMModel <- function(object, ...) {
    if (!is.element("MxRAMModel", class(object)))
    stop("\"object\" must be an object of class \"MxRAMModel\".")

    # labels of the parameters    
    my.name <- names( omxGetParameters(object) )
    # Remove NA labels
    my.name <- my.name[!is.na(my.name)]
    
    acov <- tryCatch( 2*solve(object@output$calculatedHessian)[my.name, my.name, drop=FALSE], error = function(e) e)
    # Issue a warning instead of error message
    if (inherits(acov, "error")) {
      cat("Error in solving the Hessian matrix.\n")
      warning(print(acov))
    } else {
      return(acov)
    }
}
