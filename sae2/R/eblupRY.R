     

eblupRY <- function(formula, D, T, vardir, method = c("REML","ML"),
         MAXITER = 1000, PRECISION = .1e-05, data, ...) { 
         
 if(length(method) > 1) method <- method[1]
 if (!method %in% c("REML","ML","MLE")) 
        stop(" method=\"", method, "\" must be \"REML\", or \"ML\"")
 if(class(formula)=="list") {
    NV <- length(formula)
    if(NV ==1) formula<-formula[[1]]
 } else NV <- 1
 if (class(vardir) == "list") {
   if (missing(D)) { 
     D <- length(vardir) 
   } else {
     if (D != length(vardir)) 
       stop("the length of 'vardir' must agree with D")
   }
   lapply(vardir, FUN = function(x) { if (!all(dim(x) == c(NV*T, NV*T) )) 
     stop(paste("each element of 'vardir' must be a square matrix",
       "with 'NV*T' rows"))
    } )
   vardir.temp <- matrix(0, nrow=D*NV*T, ncol=D*NV*T)
   for (d in 1:D) {
     vardir.temp[((d-1)*NV*T+1):(d*NV*T), ((d-1)*NV*T+1):(d*NV*T)] <-
       vardir[[d]]
   }
   vardir <- vardir.temp
 } else {
   if (missing(D)) stop("'D' must be specified")
   if (!is.matrix(vardir)) 
      vardir <- as.matrix(vardir)
   if (class(vardir) != "matrix") 
     stop("'vardir' must be a matrix or a list of matrices")
   if (!all(dim(vardir) == c(D*NV*T, D*NV*T) ))
     stop("'vardir' must be a square matrix with 'D*NV*T' rows")
 }
 if (NV == 1) {
    if (!missing(data)) {
        formuladata <- model.frame(formula, na.action = NULL, 
            data)
    } else {
        formuladata <- model.frame(formula, na.action = NULL)
    }
    X <- model.matrix(formula, data=formuladata)
    nformula <- nrow(X)
    if(nformula != D*T) 
       stop("length of variables must be D * T")
    y <- formuladata[,1]
    if(method=="REML") {
      result <- reml.Rao.Yu(y, X, M=D, T=T, NV=1, vcov_e=vardir,
           maxiter=MAXITER, iter.tol = PRECISION, ...) 

    } else {
      result <- mle.Rao.Yu(y, X, M=D, T=T, NV=1, vcov_e=vardir,
           maxiter=MAXITER, iter.tol = PRECISION, ...)
    }
 } else {
    depvarlist <- rep(" ", NV)
    for (nv in 1:NV) {
       formula1 <- formula[[nv]]
       depvarlist[nv] <- as.character(formula1[2])
       if (!missing(data)) {
          formuladata1 <- model.frame(formula1, na.action = NULL, 
              data)
       } else {
          formuladata1 <- model.frame(formula1, na.action = NULL)
       }
       X1 <- model.matrix(formula1, data=formuladata1)
       if (nv == 1) {
          formuladata <- list(formuladata1)
          X.list <- list(X1)
          ncol <- dim(X1)[2]
          X.names <- paste(attr(X1,"dimnames")[[2]],".1",sep="")
       } else {
          formuladata[[nv]] <- formuladata1
          X.list[[nv]] <- X1
          ncol <- ncol + dim(X1)[2]
          X.names <- append(X.names,paste(attr(X1,"dimnames")[[2]],".",nv,sep=""))
       }
    }
    y <- rep(0, D * T * NV)
    X <- matrix(0, nrow= D * T * NV, ncol = ncol)
    attr(X,"dimnames")[[2]]<-X.names
    cstart <-1
    for (nv in 1:NV) {
      X1 <- X.list[[nv]]
      y1<- formuladata[[nv]][,1]
      cend <- cstart + ncol(X1) - 1
      nformula <- nrow(X1)
      if(nformula != D*T) 
         stop("length of variables must be D * T")
      for (m in 1:D) {
        X[((m-1)*T*NV+(nv-1)*T+1):((m-1)*T*NV+nv*T),cstart:cend] <- 
           X1[((m-1)*T+1):(m*T),]
        y[((m-1)*T*NV+(nv-1)*T+1):((m-1)*T*NV+nv*T)] <-      
           y1[((m-1)*T+1):(m*T)] 
      }
      cstart <- cend + 1
    }
    if(method=="REML") {
      result <- reml.Rao.Yu(y, X, M=D, T=T, NV=NV, vcov_e=vardir,
         maxiter=MAXITER, iter.tol = PRECISION, ...) 

    } else {
      result <- mle.Rao.Yu(y, X, M=D, T=T, NV=NV, vcov_e=vardir,
           maxiter=MAXITER, iter.tol = PRECISION, ...)
    }
    colnames(result$eblup) <- depvarlist
    result$eblup <- as.data.frame(result$eblup)
    colnames(result$eblup.mse) <- depvarlist
    result$eblup.mse <- as.data.frame(result$eblup.mse)
    colnames(result$eblup.g1) <- depvarlist
    result$eblup.g1 <- as.data.frame(result$eblup.g1)
    colnames(result$eblup.g2) <- depvarlist
    result$eblup.g2 <- as.data.frame(result$eblup.g2)
    colnames(result$eblup.g3) <- depvarlist
    result$eblup.g3 <- as.data.frame(result$eblup.g3)
    colnames(result$est.fixed) <- depvarlist
    result$est.fixed <- as.data.frame(result$est.fixed)
    colnames(result$est.fixed.var) <- depvarlist
    result$est.fixed.var <- as.data.frame(result$est.fixed.var)
    colnames(result$eblup.wt1) <- depvarlist
    result$eblup.wt1 <- as.data.frame(result$eblup.wt1)
    colnames(result$eblup.wt2) <- depvarlist
    result$eblup.wt2 <- as.data.frame(result$eblup.wt2)
 }
 result$fit$iterations <- result$parm["num.iter"]
 names(result$fit$iterations) <- NULL
 std.error <- sqrt(diag(result$var.coef))
 tvalue <- result$coef/std.error
 pvalue <- ifelse (tvalue < 0, 2*pnorm(tvalue), 2*(1-pnorm(tvalue)))
 result$fit$estcoef <- data.frame(beta=result$coef, std.error=std.error,
                                  tvalue=tvalue, pvalue=pvalue,
                                  row.names=names(result$coef))
 std.error <- sqrt(diag(solve(result$inf.mat)))
 result$fit$estvarcomp <- data.frame(estimate=result$delta, 
                                     std.error=std.error,
                                     row.names=names(result$delta))
 if (method == "REML") {
   goodness <- c(result$parm["loglikelihood"], 
                 result$parm["constrained.ll"])
   names(goodness) <- c("loglike", "restrictedloglike")
 } else {
   goodness <- result$parm["loglikelihood"]
   names(goodness) <- "loglike"
 }
 result$fit$goodness <- goodness

 return(result)
}           
    
    