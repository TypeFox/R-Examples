summary.inarmix <- function(object, ...)  {
  if(object$nclasses==1)  {
    q <- length(object$coefficients) + 2
    coef.table <- matrix(0,nrow=q,ncol=4)
    coef.table[,1] <- c(object$coefficients,object$alpha,object$gamma)
    coef.table[,2] <- object$std.errs
    coef.table[,3] <- coef.table[,1]/coef.table[,2]
    coef.table[,4] <- 2*pnorm(abs(coef.table[,3]),lower.tail=F)
    coef.table[q,1] <- object$gamma + 1
    
    rownames(coef.table) <- c(object$coefnames,"autocorr.","scale")
    colnames(coef.table) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
    
    ans <- list()
    class(ans) <- "summary.inarmix"
    
    ans$call <- object$call
    ans$coef.table <- coef.table
    ans$nclasses <- 1
    ans$loglik <- object$loglikhood
    ans$bic <- object$bic
    ans$aic <- object$aic
  }
  else {
    coef.table <- list()
    mix.table <- matrix(0,nrow=object$nclasses,ncol=2)
    rownames(mix.table) <- rep("",object$nclasses)
    
    q <- ncol(object$coefficients) + 2
    for(k in 1:object$nclasses)  {
      ind1 <- (k-1)*q + 1
      ind2 <- k*q
      
      tmp <- matrix(0,nrow=q,ncol=4)
      tmp[,1] <- c(object$coefficients[k,],object$alpha[k],1 + object$gamma[k])
      tmp[,2] <- object$std.errs[ind1:ind2]
      tmp[,3] <- tmp[,1]/tmp[,2]
      tmp[q,3] <- (tmp[q,1] - 1)/tmp[q,2]
      tmp[,4] <- 2*pnorm(abs(tmp[,3]),lower.tail=F)
      
      rownames(tmp) <- c(object$coefnames,"autocorr.","scale")
      colnames(tmp) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
      coef.table[[k]] <- tmp
      
      rownames(mix.table)[k] <- paste("Group",k,sep="")
    }
    colnames(mix.table) <- c("Estimate", "Std.Err")
    mix.table[,1] <- object$mix.prop
    m1 <- q*object$nclasses + 1
    m2 <- object$nclasses*(q + 1) - 1
    mix.table[,2] <- c(object$std.errs[m1:m2],1)
    mix.table[object$nclasses,2] <- NA
    
    ans <- list()
    class(ans) <- "summary.inarmix"
    
    ans$call <- object$call
    ans$coef.table <- coef.table
    ans$mix.table <- mix.table
    ans$nclasses <- object$nclasses
    ans$niter <- object$niter
    ans$loglik <- object$loglikhood[length(object$loglikhood)]
    ans$bic <- object$bic
    ans$aic <- object$aic
    ans$reploglik <- object$reploglik
    ans$finalvals <- object$finalvals
  }  
  return(ans)
}


print.inarmix <- function(x, ...)  {
  cat("\n Call: \n",paste(deparse(x$call), sep = "\n", collapse = "\n"),"\n\n")
  
  cat("Coefficients: \n")
  print(t(x$coefficients))
  cat("\n")
  
}


print.summary.inarmix <- function(x,digits=NULL, ...)  {
  if(is.null(digits)) {
    digits <- options()$digits
  }
  else {
    options(digits = digits)
  }
  cat("\n Call: \n",paste(deparse(x$call), sep = "\n", collapse = "\n"),"\n\n")
  
  if(x$nclasses==1) {
    printCoefmat(x$coef.table, digits = digits, signif.legend=F)
    cat("\n ------------------------------ \n")          
    
    cat(apply(cbind(format(c("log-likelihood:","BIC:","AIC:"), justify="right"),
                    format(c(x$loglik,x$bic,x$aic),digits = digits), "\n"),
              1L, paste, collapse = " "), sep = "")
    cat("\n")
    invisible(x)
  }
  else {
    #### Print class proportions.
    cat("\n Class Probabilities \n")
    printCoefmat(x$mix.table, digits=digits)
    cat("\n ------------------------------ \n\n")    
    
    for(k in 1:x$nclasses)  {
      cat("Group",k," coefficients: \n\n") 
      printCoefmat(x$coef.table[[k]], digits = digits, signif.legend=F)
      
      cat("\n ------------------------------ \n\n")          
    }      
    
    
    cat(apply(cbind(format(c("log-likelihood:","BIC:","AIC:"), justify="right"),
                    format(c(x$loglik,x$bic,x$aic),digits = digits), "\n"),
              1L, paste, collapse = " "), sep = "")
    
    cat("\n Number of EM iterations: ", x$niter)
    
    cat("\n\n")
    if(length(x$reploglik) > 1)  {
      orderlogliks <- order(-x$reploglik)
      num.print <- min(5,length(orderlogliks))
      cat("Top ",num.print,"log-likelihood values: \n")
      print.obj <- matrix(x$reploglik[orderlogliks[1:num.print]],ncol=1)
      colnames(print.obj) <- c("log-lik")
      rownames(print.obj) <- rep(" ",nrow(print.obj))
      for(k in 1:nrow(print.obj))  {
        rownames(print.obj)[k] <- paste("Replicate",orderlogliks[k])
      }
      printCoefmat(print.obj,digits=digits)
      cat("\n")
    }
    invisible(x)
  }
}


diagnose <- function(results) {
  nr <- results$niter
  nclasses <- results$nclasses
  ConvMat <- matrix(0,nrow=nr,ncol=nclasses+4)
  
  ConvMat[,1] <- c(0:(nr-1))
  ConvMat[,2:(nclasses+1)] <- results$GEE.conv[1:nr,]
  ConvMat[,nclasses+2] <- results$loglikhood[1:nr]
  ConvMat[,nclasses+3] <- c(0,-diff(results$loglikhood[1:nr]))
  ConvMat[,nclasses+4] <- results$pss[1:nr]
  
  colnames(ConvMat) <- rep(" ",nclasses+4)
  colnames(ConvMat)[1] <- c("iteration")
  colnames(ConvMat)[nclasses + 2] <- c("l[i]")
  colnames(ConvMat)[nclasses + 3] <- c("l[i+1]-l[i]")
  colnames(ConvMat)[nclasses + 4] <- c("||Psi||^2")
  for(k in 2:(nclasses+1)) {
    colnames(ConvMat)[k] <- paste("GEE ",k-1,sep="") 
  }
  
  ans <- list()
  class(ans) <- "diagnose.inarmix"
  ans$em.converged <- results$em.converged
  ans$niter <- results$niter
  ans$nclasses <- results$nclasses
  ans$ConvMat <- ConvMat
  ans$loglikhood <- results$loglikhood
  return(ans)
}


print.diagnose.inarmix <- function(x,digits=NULL, ...)  {
  if(is.null(digits)) {
    digits <- options()$digits
  }
  else {
    options(digits = digits)
  }
  
  if(x$nclasses > 1)   {
    if(x$em.converged)   {
      cat("\n      Convergence was achieved. \n\n")
    }
    else {
      cat("\n      Did not converge!!!! \n\n")
    }
    
    cat("Convergence of GEE equations for
             each iteration: 1 - converged, 0 - did not converge \n\n")
    print(x$ConvMat,digits=digits)
    
    
    ll <- x$loglikhood[x$loglikhood!=0]
    tau <- length(ll)-1
    plot(0,xlim=c(0,tau),ylim=c(1.02*min(ll),.98*max(ll)),xlab="iteration",main="Log Likelihood")
    lines(c(0:tau),ll)
  }
  else   {
    stop("\n Diagnostics not applicable for one class models. \n")
  }
}