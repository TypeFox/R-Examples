###-----print.synds--------------------------------------------------------

print.synds <- function(x, ...){
  cat("Call:\n($call) ")
  print(x$call)
  cat("\nNumber of synthesised data sets: \n($m) ",x$m,"\n")  
  if (x$m==1){
    cat("\nFirst rows of synthesised data set: \n($syn)\n")
    print(head(x$syn))
  } else {
    cat("\nFirst rows of first synthesised data set: \n($syn)\n")
    print(head(x$syn[[1]]))
  }    
  cat("...\n")
  cat("\nSynthesising methods: \n($method)\n")
  print(x$method)
  cat("\nOrder of synthesis: \n($visit.sequence)\n")
  print(x$visit.sequence)
  cat("\nMatrix of predictors: \n($predictor.matrix)\n")
  print(x$predictor.matrix)     
  invisible(x)
}


###-----summary.synds------------------------------------------------------

summary.synds <- function(object, msel = NULL, 
  maxsum = 7, digits = max(3, getOption("digits")-3), ...){
  if (!is.null(msel) & !all(msel %in% (1:object$m))) stop("Invalid synthesis number(s)", call. = FALSE)

  sy <- list(m = object$m, msel = msel, method = object$method)

  if (object$m == 1){
    sy$result <- summary(object$syn,...)
  } else if (is.null(msel)){
    zall <- vector("list",object$m) 
    for (i in 1:object$m) zall[[i]] <- lapply(object$syn[[i]], summary,
      maxsum = maxsum, digits = digits, ...)
    zall.df <- Reduce(function(x,y) mapply("rbind",x,y),zall)
    meanres <- lapply(zall.df, function(x) apply(x,2,mean))
    sy$result <- summary.out(meanres)
  } else if (length(msel)==1){
    sy$result <- summary(object$syn[[msel]],...)
  } else {
    for (i in (1:length(msel))){
      sy$result[[i]] <- summary(object$syn[[msel[i]]],...)
    }
  }
  class(sy) <- "summary.synds"
  return(sy)
}


###-----print.summary.synds------------------------------------------------

print.summary.synds <- function(x, ...){

 if (x$m==1){
   cat("Synthetic object with one synethesis using methods:\n")
   print(x$method)
   cat("\n")
   print(x$result)
 } else if (is.null(x$msel)){
   cat("Synthetic object with ",x$m," syntheses using methods:\n",sep="")
   print(x$method)
   cat("\nSummary (average) for all synthetic data sets:\n",sep="")
   print(x$result)  
 } else if (length(x$msel)==1){
   cat("Synthetic object with ",x$m," syntheses using methods:\n",sep="")
   print(x$method)
   cat("\nSummary for synthetic data set ",x$msel,":\n",sep="")
   print(x$result)
 } else {
   cat("Synthetic object with ",x$m," syntheses using methods:\n",sep="")
   print(x$method)
   for (i in (1:length(x$msel))){
     cat("\nSummary for synthetic data set ",x$msel[i],":\n",sep="")
     print(x$result[[i]])
   }
 }
 invisible(x)
}


###-----mcoefvar--------------------------------------------------
# Arrange coefficients from all m syntheses in a matrix
# (same with their variances). 
# [used in lm.synds and glm.synds function]

mcoefvar <- function(analyses, ...) {
  m <- length(analyses)
  if (m == 1) {
    matcoef <- mcoefavg <- analyses[[1]]$coefficients[,1]
    matvar  <- mvaravg  <- analyses[[1]]$coefficients[,2]^2
  } else {
    namesbyfit <- lapply(lapply(analyses,coefficients),rownames)
    allnames <- Reduce(union,namesbyfit)
    matcoef <- matvar <- matrix(NA, m, length(allnames))
    dimnames(matcoef)[[2]] <- dimnames(matvar)[[2]] <- allnames
    for (i in 1:m){
      pos <- match(namesbyfit[[i]],allnames)
      matcoef[i,pos] <- analyses[[i]]$coefficients[,1]
      matvar [i,pos] <- analyses[[i]]$coefficients[,2]^2
    }
    mcoefavg <- apply(matcoef, 2, mean, na.rm = TRUE)
    mvaravg  <- apply(matvar,  2, mean, na.rm = TRUE)
    #bm <- apply(matcoef,2,var) not needed xpt for partial synthesis
  }
  if (m > 1) rownames(matcoef) <- rownames(matvar) <- paste0("syn=", 1:m)
  return(list(mcoef    = matcoef,  mvar    = matvar, 
              mcoefavg = mcoefavg, mvaravg = mvaravg))
}


###-----lm.synds-----------------------------------------------------------

lm.synds <- function(formula, data, ...)
{
 if (!class(data)=="synds") stop("Data must have class synds\n")
 call <- match.call()
 fitting.function <- "lm"
 analyses <- as.list(1:data$m)

 # do the repated analysis, store the result without data
 if (data$m==1) {
   analyses[[1]] <- summary(lm(formula,data=data$syn,...))
 } else {
   for(i in 1:data$m) {
     analyses[[i]] <- summary(lm(formula,data=data$syn[[i]],...))
   }
 }
 allcoefvar <- mcoefvar(analyses = analyses)
      
 # return the complete data analyses as a list of length m
 object <- list(call=call, mcoefavg = allcoefvar$mcoefavg, 
             mvaravg = allcoefvar$mvaravg, proper = data$proper, m = data$m, 
             analyses = analyses, fitting.function = fitting.function,
             n = data$n, k=data$k, mcoef = allcoefvar$mcoef,
             mvar = allcoefvar$mvar)
 class(object) <- "fit.synds"
 return(object)
}


###-----glm.synds----------------------------------------------------------

glm.synds <- function(formula, family="binomial", data, ...)
{
 if (!class(data)=="synds") stop("Data must have class synds\n")
 call <- match.call()
 fitting.function <- "glm"
 analyses <- as.list(1:data$m)
 
 # do the repated analysis, store the result without data
 if (data$m==1) {
   analyses[[1]] <- summary(glm(formula,data=data$syn,family=family,...))
 } else {
   for(i in 1:data$m) {
     analyses[[i]] <- summary(glm(formula,data=data$syn[[i]],family=family,...))
   }
 }
 allcoefvar <- mcoefvar(analyses = analyses)
 
 # return the complete data analyses as a list of length m
 object <- list(call=call, mcoefavg = allcoefvar$mcoefavg, 
             mvaravg = allcoefvar$mvaravg, proper = data$proper, m = data$m,
             analyses = analyses, fitting.function = fitting.function,
             n = data$n, k = data$k, mcoef = allcoefvar$mcoef,
             mvar = allcoefvar$mvar)
 class(object) <- "fit.synds"
 return(object)
}


###-----print.fit.synds----------------------------------------------------

print.fit.synds <- function(x, msel = NULL, ...)
{
  if (!is.null(msel) & !all(msel %in% (1:x$m))) stop("Invalid synthesis number(s)", call. = FALSE)
  cat("\nCall:\n")
  print(x$call)
  if (is.null(msel)){
    cat("\nCombined coefficient estimates:\n")
    print(x$mcoefavg)
  } else {
    cat("\nCoefficient estimates for selected syntheses:\n")
    print(x$mcoef[msel,,drop=FALSE])
  }
  invisible(x)
}


###-----summary.fit.synds--------------------------------------------------

summary.fit.synds <- function(object, population.inference = FALSE, msel = NULL, partly = FALSE, ...)
{ # df.residual changed to df[2] because didn't work for lm - check if that's ok
  if (!class(object) == "fit.synds") stop("Object must have class fit.synds\n")
  m <- object$m
  k <- object$k
  n <- object$n
  
  coefficients <- object$mcoefavg     # mean of coefficients (over m syntheses)
  vars <- object$mvaravg              # mean of variances (over m syntheses)

  if (population.inference == F){ ## inf to Q hat

    if(partly == TRUE){
      bm = c(0)
      for(i in 1:m){
        bm = bm + ((object$mcoef[i,] - object$mcoefavg)^2 / (m - 1))
      }
      result <- cbind(coefficients,
                      sqrt(bm/m),
                      sqrt(vars*k/n),
                      coefficients/sqrt(vars*k/n),
                      #sqrt((1 + coefficients^2/vars/2/object$analyses[[1]]$df[2])*n/k/m)
                      sqrt(1 + coefficients^2/vars/4/object$analyses[[1]]$df[2] * n/k/m))  #!?????
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","se(Beta).syn","Z.syn","se(Z.syn)")
    }
    else if (object$proper == F){
      ## simple synthesis
      result <- cbind(coefficients,
                      sqrt(vars/m),
                      sqrt(vars*k/n),
                      coefficients/sqrt(vars*k/n),
                      #sqrt((1 + coefficients^2/vars/2/object$analyses[[1]]$df[2])*n/k/m)
                      sqrt(1 + coefficients^2/vars/4/object$analyses[[1]]$df[2] * n/k/m))  #!?????
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","se(Beta).syn","Z.syn","se(Z.syn)")
    } else {
      ## proper synthesis
      result <- cbind(coefficients,
                      sqrt(vars*(1+k/n)/m), 
                      sqrt(vars*k/n),
                      coefficients/sqrt(vars*k/n),
                      #sqrt((1 + k/n + coefficients^2/vars/2/object$analyses[[1]]$df[2])*n/k/m)
                      sqrt(1 + coefficients^2/vars/4/object$analyses[[1]]$df[2]*n/k/m))    #!?????
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","se(Beta).syn","Z.syn","se(Z.syn)")
    }

  } else { ## pop inference to Q

  	if(partly == TRUE){
      bm = c(0)
      for(i in 1:m){
        bm = bm + ((object$mcoef[i,] - object$mcoefavg)^2 / (m - 1))
      }
      result <- cbind(coefficients,
                      sqrt(bm/m + vars),
                      coefficients/sqrt(bm/m + vars))
                      #sqrt((1 + coefficients^2/vars/2/object$analyses[[1]]$df[2])*n/k/m)
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","Z.syn")
    }
    else if (object$proper == F){
      ## simple synthesis
      Tf <- vars*(k/n+1/m)
      result <- cbind(coefficients,sqrt(Tf),coefficients/sqrt(Tf))
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","Z.syn")
    }
    else {
      ## proper synthesis
      Tf <- vars*(k/n+(1+k/n)/m)
      result <- cbind(coefficients,sqrt(Tf),coefficients/sqrt(Tf))
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","Z.syn")
    }
  }
  res <- list(call = object$call, proper = object$proper,
              population.inference = population.inference,
              fitting.function = object$fitting.function,
              m = m, coefficients = result, n = n, k = k, 
              analyses = object$analyses, msel = msel)
  class(res) <- "summary.fit.synds"
  return(res)
}


###-----print.summary.fit.synds--------------------------------------------

print.summary.fit.synds <- function(x, ...) {
  if (!is.null(x$msel) & !all(x$msel %in% (1:x$m))) stop("Invalid synthesis number(s)", call. = FALSE)
  
  if (is.null(x$msel)){
    if (x$m==1) {
      cat("\nFit to synthetic data set with a single synthesis.\n")
    } else {
      cat("\nFit to synthetic data set with ",x$m," syntheses.\n",sep="")
    }
    if (x$population.inference) {
      cat("Inference to population coefficients.\n")
    } else {
      cat("Inference to coefficients and standard errors\nthat would be obtained from the observed data.\n")
    }
    cat("\nCall:\n")
    print(x$call)
    cat("\nCombined estimates:\n")
    if (x$population.inference){
      print(x$coefficients[,c("B.syn","se(B.syn)","Z.syn")])
    } else {
      print(x$coefficients[,c("B.syn","se(Beta).syn","Z.syn")])
    }
  } else {
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficient estimates for selected syntheses:\n")
    for(i in x$msel) {          
      cat("\nsyn=",i,"\n",sep="")
      print(x$analyses[[i]]$coefficients)
    }
  }      
  invisible(x)
}


###-----print.compare.fit.synds--------------------------------------------

print.compare.fit.synds <- function(x, ...){
  cat("\nCall used to fit models to the data:\n")
  print(x$call)
  if (!is.null(x$coef.obs)){
    cat("\nEstimates for the observed data set:\n")
    print(x$coef.obs)
    cat("\nCombined estimates for the synthetised data set(s):\n")
    print(x$coef.syn[,c("B.syn","se(Beta).syn","se(B.syn)","Z.syn")])
    cat("\nDifferences synthetic vs. observed:\n")
    print(cbind.data.frame(x$coef.diff,x$ci.overlap))                                 
  }
  if (!is.null(x$ci.plot)){
    cat("\nConfidence interval plot:\n")
    print(x$ci.plot)
  }
  invisible(x)
}


###-----print.compare.synds------------------------------------------------

print.compare.synds <- function(x, ...) {
  cat("\nComparing percentages observed with synthetic\n\n")
  if (class(x$plots)[1]=="gg"){
    print(x$tables) 
    print(x$plots)
  } else {
    for (i in 1:length(x$tables)) {
      print(x$tables[[i]]) 
      print(x$plots[[i]])
      if (i < length(x$tables)) {
        cat("Press return for next plot: ")
        ans <- readline()
      }
    }
  }
 invisible(x)
}


###-----print.utility.synds------------------------------------------------

print.utility.synds <- function(x, ...){
  
	cat("\nMethod:\n($method)\n")
	cat(x$method,"\n")
	
	cat("\nMean and SD of propensity score utility value: \n($utility.summary)\n")
	print(x$utility.summary)
	
#	cat("\nPropensity score utility values for each synthetic set \n(first few if numerous): \n($utility.raw)\n")
#	print(head(x$utility.raw))

  if(!is.null(x$deviance.summary)){
    cat("\nMean and SD of propensity deviance statistics: \n($deviance.summary)\n")
    print(x$deviance.summary)

#    cat("\nPropensity deviance statistics for each synthetic set \n(first few if numerous): \n($deviance.raw)\n")
#    print(head(x$deviance.raw))
  }
	
	if(!is.null(x$nullstats.summary)){
    cat("\nNull distribution statistics: \n($nullstats.summary)\n")
    print(x$nullstats.summary)
	}
	
	invisible(x)
}

                               

###-----summary.out--------------------------------------------------------
summary.out <- function (z, digits = max(3L, getOption("digits") - 3L), ...)
{
    ncw <- function(x) {
        zz <- nchar(x, type = "w")
        if (any(na <- is.na(zz))) {
            zz[na] <- nchar(encodeString(zz[na]), "b")
        }
        zz
    }
    nv <- length(z)
    nm <- names(z)
    lw <- numeric(nv)
    nr <- if (nv)
        max(unlist(lapply(z, NROW)))
    else 0
    for (i in seq_len(nv)) {
        sms <- z[[i]]
        if (is.matrix(sms)) {
            cn <- paste(nm[i], gsub("^ +", "", colnames(sms),
                useBytes = TRUE), sep = ".")
            tmp <- format(sms)
            if (nrow(sms) < nr)
                tmp <- rbind(tmp, matrix("", nr - nrow(sms),
                  ncol(sms)))
            sms <- apply(tmp, 1L, function(x) paste(x, collapse = "  "))
            wid <- sapply(tmp[1L, ], nchar, type = "w")
            blanks <- paste(character(max(wid)), collapse = " ")
            wcn <- ncw(cn)
            pad0 <- floor((wid - wcn)/2)
            pad1 <- wid - wcn - pad0
            cn <- paste0(substring(blanks, 1L, pad0), cn, substring(blanks,
                1L, pad1))
            nm[i] <- paste(cn, collapse = "  ")
            z[[i]] <- sms
        }
        else {
            sms <- format(sms, digits = digits)
            lbs <- format(names(sms))
            sms <- paste0(lbs, ":", sms, "  ")
            lw[i] <- ncw(lbs[1L])
            length(sms) <- nr
            z[[i]] <- sms
        }
    }
    if (nv) {
        z <- unlist(z, use.names = TRUE)
        dim(z) <- c(nr, nv)
        if (anyNA(lw))
            warning("probably wrong encoding in names(.) of column ",
                paste(which(is.na(lw)), collapse = ", "))
        blanks <- paste(character(max(lw, na.rm = TRUE) + 2L),
            collapse = " ")
        pad <- floor(lw - ncw(nm)/2)
        nm <- paste0(substring(blanks, 1, pad), nm)
        dimnames(z) <- list(rep.int("", nr), nm)
    }
    else {
        z <- character()
        dim(z) <- c(nr, nv)
    }
    attr(z, "class") <- c("table")
    z
}

