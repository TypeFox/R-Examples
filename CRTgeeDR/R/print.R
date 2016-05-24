### summary function for CRTgeeDR object.
#' Summarizing CRTgeeDR object.
#' 
#' Summary CRTgeeDR object
#' 
#' @param object CRTgeeDR object
#' @param ... ignored
#' 
#' @aliases summary.CRTgeeDR summary
#' @method summary CRTgeeDR
#' @export
summary.CRTgeeDR <- function(object, ...)  {
  Coefs <- matrix(0,nrow=length(object$beta),ncol=5)
  Coefs[,1] <- c(object$beta)
  naive <- is.character(object$var)
  Coefs[,2] <- sqrt(diag(object$var.naiv))
  if(naive){Coefs[,3] <- rep(0, length(object$beta))}else{Coefs[,3] <- sqrt(diag(object$var))}
  if(naive){Coefs[,4] <- Coefs[,1]/Coefs[,2]}else{Coefs[,4] <- Coefs[,1]/Coefs[,3]}
  Coefs[,5] <- round(2*pnorm(abs(Coefs[,4]), lower.tail=F), digits=8)
  colnames(Coefs) <- c("Estimates","Naive SE","Robust SE", "wald", "p")
  
  summ <- list(beta = Coefs[,1], se.model = Coefs[,2], se.robust = Coefs[,3], wald.test = Coefs[,4], p = Coefs[,5],
               alpha = object$alpha, corr = object$corr, phi = object$phi, niter = object$niter, clusz = object$clusz, coefnames = object$coefnames, weights=object$weights)
  
  class(summ) <- 'summary.CRTgeeDR'
  return(summ)
}

### summary function for CRTgeeDR object.
#' Print the summarizing CRTgeeDR object.
#' 
#' Print Summary CRTgeeDR object
#' 
#' @param x summary.CRTgeeDR x
#' @param ... ignored
#' 
#' @aliases print.summary.CRTgeeDR print.summary
#' @method print summary.CRTgeeDR
#' @export
print.summary.CRTgeeDR <- function(x, ...){
  Coefs <- matrix(0,nrow=length(x$coefnames),ncol=5)
  rownames(Coefs) <- c(x$coefnames)
  colnames(Coefs) <- c("Estimates","Model SE","Robust SE", "wald", "p")
  Coefs[,1] <- x$beta
  Coefs[,2] <- x$se.model
  Coefs[,3] <- x$se.robust
  Coefs[,4] <- x$wald.test
  Coefs[,5] <- x$p
  
  #print("Call: ", object$call, "\n")
  print(signif(Coefs, digits=4))
  cat("\n Est. Correlation: ", signif(x$alpha, digits=4), "\n")
  cat(" Correlation Structure: ", x$corr, "\n")
  cat(" Est. Scale Parameter: ", signif(x$phi, digits=4), "\n")
  cat("\n Number of GEE iterations:", x$niter, "\n")
  cat(" Number of Clusters: ", length(x$clusz), "   Maximum Cluster Size: ", max(x$clusz), "\n")
  cat(" Number of observations with nonzero weight: ", sum(x$weights != 0), "\n")
}



#' Prints CRTgeeDR object.
#' 
#' Prints CRTgeeDR object
#' 
#' @param x CRTgeeDR x
#' @param ... ignored
#' 
#' @aliases print.CRTgeeDR print
#' @method print CRTgeeDR
#' @export
### print function for CRTgeeDR object
print.CRTgeeDR <- function(x, ...){
  coefdf <- signif(data.frame(x$beta), digits=4)
  rownames(coefdf) <- x$coefnames
  colnames(coefdf) <- ""
  print(x$call)
  cat("\n", "Coefficients:", "\n")
  print(t(coefdf))
  cat("\n Scale Parameter: ", signif(x$phi, digits=4), "\n")
  cat("\n Correlation Model: ", x$corr)
  cat("\n Estimated Correlation Parameters: ", signif(x$alpha, digits=4), "\n")
  cat("\n Number of clusters: ", length(x$clusz), "  Maximum cluster size: ", max(x$clusz), "\n")
  cat(" Number of observations with nonzero weight: ", sum(x$weights != 0), "\n")
}

#' Get Mean, Sd and CI for estimates from CRTgeeDR object.
#' 
#' Get the estimates, standard deviations and confidence intervals from an CRTgeeDR object associated with a regressor given in argument.
#' 
#' @param object CRTgeeDR
#' @param nameTRT, character including the name of the variable of interest (often the treatment)
#' @param quantile, value of the normal quantile for the IC. default is 1.96 for 95\%CI. 
#' @export 
#' 
getCI <-function(object,nameTRT="TRT",quantile=1.96){
  stats.summary <- c()
  stats.summary<-c(stats.summary,object$beta[which(object$coefnames==nameTRT)])
 
  stats.summary<-c(stats.summary,sqrt(object$var.naiv[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)]))
  stats.summary<-c(stats.summary,object$beta[which(object$coefnames==nameTRT)]-quantile*sqrt(object$var.naiv[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)]))
  stats.summary<-c(stats.summary,object$beta[which(object$coefnames==nameTRT)]+quantile*sqrt(object$var.naiv[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)]))
  
  stats.summary<-c(stats.summary,ifelse(is.null(object$var),NA,sqrt(object$var[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)])))
  stats.summary<-c(stats.summary,ifelse(is.null(object$var),NA,object$beta[which(object$coefnames==nameTRT)]-quantile*sqrt(object$var[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)])))
  stats.summary<-c(stats.summary,ifelse(is.null(object$var),NA,object$beta[which(object$coefnames==nameTRT)]+quantile*sqrt(object$var[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)])))
  
  stats.summary<-c(stats.summary,ifelse(is.null(object$var.nuisance),NA,sqrt(object$var.nuisance[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)])))
  stats.summary<-c(stats.summary,ifelse(is.null(object$var.nuisance),NA,object$beta[which(object$coefnames==nameTRT)]-quantile*sqrt(object$var.nuisance[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)])))
  stats.summary<-c(stats.summary,ifelse(is.null(object$var.nuisance),NA,object$beta[which(object$coefnames==nameTRT)]+quantile*sqrt(object$var.nuisance[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)])))

  stats.summary<-c(stats.summary,ifelse(is.null(object$var.fay),NA,sqrt(object$var.fay[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)])))
  stats.summary<-c(stats.summary,ifelse(is.null(object$var.fay),NA,object$beta[which(object$coefnames==nameTRT)]-quantile*sqrt(object$var.fay[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)])))
  stats.summary<-c(stats.summary,ifelse(is.null(object$var.fay),NA,object$beta[which(object$coefnames==nameTRT)]+quantile*sqrt(object$var.fay[which(object$coefnames==nameTRT),which(object$coefnames==nameTRT)])))
  
  stats.summary<-c(stats.summary,object$converged)
  stats.summary<-t(as.matrix(stats.summary))
  colnames(stats.summary)<-c("Estimate","Naive SD","CI naive min","CI naive max","Sandwich SD","CI Sandwich min","CI Sandwich max","Nuisance-adj SD","CI Nuisance-adj min","CI Nuisance-adj max","Fay-adj SD","CI Fay-adj min","CI Fay-adj max","Convergence Status")
  rownames(stats.summary)<-c(nameTRT)
  return(stats.summary)
}

#' Get the histogram of weights for IPW and adequation for the glm weights model
#' 
#' Get the histogram and some basic statistics for the weights used in the IPW part.
#' 
#' @param object CRTgeeDR
#' @param save, logical if TRUE the plot is saved as a pdf in the current directory
#' @param name, name of the plot saved as pdf
#' @param typeplot, integer indicating which is the adequation diagnostic plot for the PS. Default is NULL no output. '0', all available in plot.glm are displayed, '1' Residuals vs Fitted, '2' Normal Q-Q,
#'        '3' Scale-Location, '4' Cook's distance, '5' Residuals vs Leverage and '6' Cook's dist vs Leverage* h[ii] / (1 - h[ii])
#' @export 
#' 
getPSPlot <-function(object,save=FALSE,name="plotPS",typeplot=NULL){
  ..count.. <- NULL
  if(save){pdf(paste(name,".pdf",sep=""))}
    toplot<-as.data.frame(object$weights)
    names(toplot)<-"weights"
    m <- ggplot(toplot, aes(x=weights))
    plot<-m + geom_histogram(aes(y=..count..,fill=..count..),binwidth=0.1)+
      ggtitle(paste("Summary: Q1=",round(quantile(object$used.weights,0.25),2)," Q2=",round(quantile(object$used.weights,0.5),2)," Q3=",round(quantile(object$used.weights,0.75),2)," max=",round(max(object$used.weights),2), sep=' '))  
    plot(plot)
    if(!is.null(object$ps.model)){
    if(!is.null(typeplot)){
     if(typeplot==0){
       for (i in 1:6) {
         plot(object$ps.model,i)
       }
     }else{
       plot(object$ps.model,typeplot)
     }
   }
  }else{
    stop("Adequation plots not available. PS had not been computed internally")
  }
  if(save){dev.off()}
    toplot<-as.data.frame(object$weights)
    names(toplot)<-"weights"
    m <- ggplot(toplot, aes(x=weights))
    plot<-m + geom_histogram(aes(y=..count..,fill=..count..),binwidth=0.1)+
      ggtitle(paste("Summary: Q1=",round(quantile(object$used.weights,0.25),2)," Q2=",round(quantile(object$used.weights,0.5),2)," Q3=",round(quantile(object$used.weights,0.75),2)," max=",round(max(object$used.weights),2), sep=' '))  
    plot(plot)
    if(!is.null(object$ps.model)){
      if(!is.null(typeplot)){
        if(typeplot==0){
          for (i in 1:6) {
            plot(object$ps.model,i)
          }
        }else{
          plot(object$ps.model,typeplot)
        }
      }
    }else{
      stop("Adequation plots not available. PS had not been computed internally")
    }
}



#' Get the observed vs fitted residuals
#' 
#' Get the histogram and some basic statistics for the weights used in the IPW part.
#' 
#' @param object CRTgeeDR
#' @param save, logical if TRUE the plot is saved as a pdf in the current directory
#' @param name, name of the plot saved as pdf
#' @param typeplot, integer indicating which is the adequation diagnostic plot for the PS. '0', all available in plot.glm are displayed, '1' Residuals vs Fitted, '2' Normal Q-Q,
#'        '3' Scale-Location, '4' Cook's distance, '5' Residuals vs Leverage and '6' Cook's dist vs Leverage* h[ii] / (1 - h[ii])
#' @export 
#' 
getOMPlot <-function(object,save=FALSE,name="plotOM",typeplot=0){
    if(is.null(object$call$model.augmentation.trt)&(is.null(object$call$aug))){
      stop("The object given as attributes is not the result of an AUGMENTED analysis")
    } else{
      if(!is.null(object$call$model.augmentation.trt)){
          if(typeplot==0){
          if(save){pdf(paste(name,"_trt.pdf",sep=""))}
          print("For the treated group:")
          plot(object$om.model.trt,1:6)
          if(save){dev.off()}
          
          if(save){pdf(paste(name,"_ctrl.pdf",sep=""))}
          print("For the control group:")
          plot(object$om.model.ctrl,1:6)
          if(save){dev.off()}
        }else{
          print("For the treated group:")
          if(save){pdf(paste(name,"_trt.pdf",sep=""))}
          plot(object$om.model.trt,typeplot)
          if(save){dev.off()}
          
          if(save){pdf(paste(name,"_ctrl.pdf",sep=""))}
          print("For the control group:")
          plot(object$om.model.ctrl,typeplot)
          if(save){dev.off()}
          
        }
      }else{
        stop("Adequation plots not available. OM had not been computed internally")
      }
    }

}



