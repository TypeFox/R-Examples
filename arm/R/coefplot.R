coefplot.default <- function(coefs, sds, CI=2, 
            lower.conf.bounds, upper.conf.bounds,
            varnames=NULL, 
            vertical=TRUE,
            v.axis=TRUE, h.axis=TRUE,
            cex.var=0.8, cex.pts=0.9, col.pts=1, pch.pts=20,
            var.las=2, main=NULL, xlab=NULL, ylab=NULL, mar=c(1,3,5.1,2),
            plot=TRUE, add=FALSE, offset=0.1, ...)
{
  
     # collect informations
    if (is.list(coefs)){
      coefs <- unlist(coefs)
    }
    
    
    n.x <- length(coefs)
    idx <- seq(1, n.x)   
    #bound <- lower.bound
    if(!missing(lower.conf.bounds)){
      if(length(coefs)!=length(lower.conf.bounds)){
        stop("Number of conf.bounds does not equal to number of estimates")
      }
    }
    if(!missing(upper.conf.bounds)){
      if(length(coefs)!=length(upper.conf.bounds)){
        stop("Number of conf.bounds does not equal to number of estimates")
      }
    }

    if(!missing(sds)){
      coefs.h <- coefs + CI*sds 
      coefs.l <- coefs - CI*sds
      est1 <- cbind(coefs - sds, coefs + sds)
      est2 <- cbind(coefs - 2*sds, coefs + 2*sds)
      if(!missing(lower.conf.bounds)){
        est1[,1] <- lower.conf.bounds
        CI <- 1
      }
      if(!missing(upper.conf.bounds)){
        est1[,2] <- upper.conf.bounds
        CI <- 1
      }

    }else{
      #coefs.h <- upper.conf.bounds
      #coefs.l <- lower.conf.bounds
      est1 <- cbind(coefs, coefs)
      if(!missing(lower.conf.bounds)){
        est1[,1] <- lower.conf.bounds
        CI <- 1
      }
      if(!missing(upper.conf.bounds)){
        est1[,2] <- upper.conf.bounds
        CI <- 1
      }
    }
    old.par <- par(no.readonly=TRUE)
    #on.exit(par(old.par))  
    min.mar <- par('mar')
    
    if (is.null(main)){main <- "Regression Estimates"}
    if (is.null(xlab)){xlab <- ""}
    if (is.null(ylab)){ylab <- ""}
        
    par(mar = mar)
    
    if (is.null(varnames)) {
      maxchar <- 0
    }
    else{
      maxchar <- max(sapply(varnames, nchar))
    }

    
    # add margin to the axis
    k <- 1/n.x   
    if(plot){
      if (vertical){

        mar[2] <- max(min.mar[2], trunc(mar[2] + maxchar/10)) + 0.1

        par(mar=mar)
        if(!add){
          plot(c(coefs.l, coefs.h), c(idx+k,idx-k), type="n",                                     
            axes=F, main=main, xlab=xlab, ylab=ylab,...) 
          if (h.axis){                                                  
            #axis(1)                                
            axis(3)
          }
          if (v.axis){
            axis(2, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
              lty=0, cex.axis=cex.var) 
          }
          abline(v=0, lty=2)                                                 
          points(coefs, idx, pch=pch.pts, cex=cex.pts, col=col.pts)
          if (CI==2){
            segments (est1[,1], idx, est1[,2], idx, lwd=2, col=col.pts)     
            segments (est2[,1], idx, est2[,2], idx, lwd=1, col=col.pts)
          }
          else{
            segments (est1[,1], idx, est1[,2], idx, lwd=1, col=col.pts)     
          }
        }
        else{
          idx <- idx + offset
          points(coefs, idx, pch=pch.pts, cex=cex.pts, col=col.pts)
          if (CI==2){
            segments (est1[,1], idx, est1[,2], idx, lwd=2, col=col.pts)     
            segments (est2[,1], idx, est2[,2], idx, lwd=1, col=col.pts)
          }
          else{
            segments (est1[,1], idx, est1[,2], idx, lwd=1, col=col.pts)     
          }
        }
    } # end of if vertical
    else{ # horizontal
      mar[1] <- max(min.mar[1], trunc(mar[1] + maxchar/10)) + 0.1
      par(mar=mar)
      if(!add){
        plot(c(idx+k,idx-k), c(coefs.l, coefs.h), type="n", axes=F, 
          main=main, xlab=xlab, ylab=ylab,...)                                                  
        if (v.axis){
          axis(2, las=var.las)                                
          #axis(4, las=var.las)
        }
        if (h.axis){
          axis(1, 1:n.x, varnames[1:n.x], las=var.las, tck=FALSE, 
            lty=0, cex.axis=cex.var) 
        }
        abline(h=0, lty=2)                                                 
        points(idx, coefs, pch=pch.pts, cex=cex.pts, col=col.pts)
        if (CI==2){
          segments (idx, est1[,1], idx, est1[,2], lwd=2, col=col.pts)     
          segments (idx, est2[,1], idx, est2[,2], lwd=1, col=col.pts)
        }
        else if (CI==1) {
          segments (idx, est1[,1], idx, est1[,2], lwd=1, col=col.pts)     
        }
      }
      else{
        idx <- idx + offset
        points(idx, coefs, pch=pch.pts, cex=cex.pts, col=col.pts)
        if (CI==2){
          segments (idx, est1[,1], idx, est1[,2], lwd=2, col=col.pts)     
          segments (idx, est2[,1], idx, est2[,2], lwd=1, col=col.pts)
        }
        else if (CI==1) {
          segments (idx, est1[,1], idx, est1[,2], lwd=1, col=col.pts)     
        }
      }
    }   
  }
  else{
    if (vertical){
      mar[2] <- max(min.mar[2], trunc(mar[2] + maxchar/10))  + 0.1
      par(mar=mar)
      plot(c(coefs.l, coefs.h), c(idx+k,idx-k), type="n",                                     
          axes=F, main="", xlab=xlab, ylab=ylab,...)
#      if (v.axis){
#          axis(2, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
#              lty=0, cex.axis=cex.var) 
#      }
    }
    else{ # horizontal
      mar[1] <- max(min.mar[1], trunc(mar[1] + maxchar/10)) + 0.1
      par(mar=mar)
      plot(c(idx+k,idx-k), c(coefs.l, coefs.h), type="n", axes=F, 
        main=main, xlab=xlab, ylab=ylab,...)                                                  
      #if (h.axis){
#          axis(1, 1:n.x, varnames[1:n.x], las=var.las, tck=FALSE, 
#              lty=0, cex.axis=cex.var) 
#      }
    }
  } 
  #on.exit(par(old.par))  
}

setMethod("coefplot", signature(object = "numeric"),
  function(object, ...)
{
  coefplot.default(object, ...)
} 
)



setMethod("coefplot", signature(object = "lm"), 
    function(object, varnames=NULL, intercept=FALSE, ...)
    {
    # collect informations
    coefs <- summary(object)$coef[,1]
    sds <- summary(object)$coef[,2]
    ifelse (is.null(varnames), varnames <- names(coefs),
            varnames <- varnames)
    if (length(varnames)!= length(names(coefs))){
      stop(message="the length of varnames does not equal the length of predictors.  
      Note: varnames must include a name for constant/intercept")
    }
    chk.int <- attr(terms(object), "intercep")
    if(chk.int & intercept | !chk.int & intercept | !chk.int & !intercept){
      intercept <- TRUE
      coefs <- coefs
      sds <- sds
      varnames <- varnames
    } else if(chk.int & !intercept){
      coefs <- coefs[-1]
      sds <- sds[-1]
      varnames <- varnames[-1]
    }    
    
    
    # plotting
    coefplot(coefs, sds, 
        varnames=varnames, ...)
    }
)
           
setMethod("coefplot", signature(object = "glm"),
    function(object, varnames=NULL, intercept=FALSE,...)
    {
    # collect informations
    coefs <- summary(object)$coef[,1]
    sds <- summary(object)$coef[,2]
    ifelse (is.null(varnames), varnames <- names(coefs),
            varnames <- varnames)
    if (length(varnames)!= length(names(coefs))){
      stop(message="the length of varnames does not equal the length of predictors.  
      Note: varnames must include a name for constant/intercept")
    }    
    chk.int <- attr(terms(object), "intercep")
    if(chk.int & intercept | !chk.int & intercept | !chk.int & !intercept){
      intercept <- TRUE
      coefs <- coefs
      sds <- sds
      varnames <- varnames
    } else if(chk.int & !intercept){
      coefs <- coefs[-1]
      sds <- sds[-1]
      varnames <- varnames[-1]
    }    
    
    # plotting
    coefplot(coefs, sds, 
        varnames=varnames, ...)
    }                                                                         
)


setMethod("coefplot", signature(object = "bugs"),
    function(object, var.idx=NULL, varnames=NULL, 
            CI=1, vertical=TRUE,
            v.axis=TRUE, h.axis=TRUE, 
            cex.var=0.8, cex.pts=0.9, 
            col.pts=1, pch.pts=20, var.las=2, 
            main=NULL, xlab=NULL, ylab=NULL, 
            plot=TRUE, add=FALSE, offset=.1,
            mar=c(1,3,5.1,2), ...)
{  
    
  if (is.null(var.idx)){
    var.idx <- 1:length(object$summary[,"50%"])
  }
  n.x <- length(var.idx)
  idx <- 1:n.x
  
  coefs <- object$summary[,"50%"][var.idx]
  if (is.null(varnames)){
    varnames <- names(coefs)     
  }
  
  if (is.null(main)){main <- "Regression Estimates"}
  if (is.null(xlab)){xlab <- ""}
  if (is.null(ylab)){ylab <- ""}
  
  min.mar <- par('mar')  
  par(mar=mar)
  
  
  maxchar <- max(sapply(varnames, nchar))

  
  k <- 1/n.x
  
  if (CI==1){
    CI50.h <- object$summary[,"75%"][var.idx]
    CI50.l <- object$summary[,"25%"][var.idx]
    CI50 <- cbind(CI50.l, CI50.h)
    if (vertical){
      mar[2] <- min(min.mar[2], trunc(mar[2] + maxchar/10)) + 0.1
      par(mar=mar)
      if(add){
        segments (CI50[,1], idx+offset, CI50[,2], idx+offset, lwd=1, col=col.pts)     
        points(coefs, idx+offset, pch=20, cex=cex.pts, col=col.pts)
      }
      else{
        plot(c(CI50[,1],CI50[,2]), c(idx+k,idx-k), type="n", 
          axes=F, main=main, xlab=xlab, ylab=ylab, ...) 
        if(plot){
          if (h.axis){
            axis(3)
          }
          if (v.axis){
          axis(2, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
              lty=0, cex.axis=cex.var)  
          }
          abline(v=0, lty=2)                                                 
          segments (CI50[,1], idx, CI50[,2], idx, lwd=1, col=col.pts)     
          points(coefs, idx, pch=20, cex=cex.pts, col=col.pts)  
        }
      }
    }
    else {
      mar[1] <- min(min.mar[1], trunc(mar[1] + maxchar/10))  + 0.1
      par(mar=mar)
      if(add){
          segments (idx+offset, CI50[,1], idx+offset, CI50[,2], lwd=1, col=col.pts)     
          points(idx+offset, coefs, pch=20, cex=cex.pts, col=col.pts)
      }
      else{
        plot(c(idx+k,idx-k), c(CI50[,1],CI50[,2]), type="n",                                     
            axes=F, main=main, xlab=xlab, ylab=ylab,...) 
        if(plot){
          if (v.axis){
            axis(2)
          }
          if (h.axis){
            axis(1, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
              lty=0, cex.axis=cex.var)  
          }
          
          abline(h=0, lty=2)                                                 
          segments (idx, CI50[,1], idx, CI50[,2], lwd=1, col=col.pts)     
          points(idx, coefs, pch=20, cex=cex.pts, col=col.pts)
        }
      }
    }   
  }
  
  if (CI==2){
    CI50.h <- object$summary[,"75%"][var.idx]
    CI50.l <- object$summary[,"25%"][var.idx]
    CI95.h <- object$summary[,"97.5%"][var.idx]
    CI95.l <- object$summary[,"2.5%"][var.idx]
    CI50 <- cbind(CI50.l, CI50.h)
    CI95 <- cbind(CI95.l, CI95.h)
    if (vertical){
      mar[2] <- min(min.mar[2], trunc(mar[2] + maxchar/10)) + 0.1
      par(mar=mar)
      if(add){
        segments (CI50[,1], idx+offset, CI50[,2], idx+offset, lwd=2, col=col.pts) 
        segments (CI95[,1], idx+offset, CI95[,2], idx+offset, lwd=1, col=col.pts)    
        points(coefs, idx+offset, pch=20, cex=cex.pts, col=col.pts)
      }
      else{
        plot(c(CI95[,1],CI95[,2]), c(idx+k,idx-k), type="n",                                     
          axes=F, main=main, xlab=xlab, ylab=ylab,...) 
        if(plot){
          if (h.axis){
            axis(3)
          }
          if (v.axis){
            axis(2, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
              lty=0, cex.axis=cex.var)  
          }
          abline(v=0, lty=2)                                                 
          segments (CI50[,1], idx, CI50[,2], idx, lwd=2, col=col.pts) 
          segments (CI95[,1], idx, CI95[,2], idx, lwd=1, col=col.pts)    
          points(coefs, idx, pch=20, cex=cex.pts, col=col.pts)
        }
      }
    }
    else {
      mar[1] <- min(min.mar[1], trunc(mar[1] + maxchar/10)) + 0.1
      par(mar=mar)
      if(add){
        segments (idx+offset, CI50[,1], idx+offset, CI50[,2], lwd=2, col=col.pts)
        segments (idx+offset, CI95[,1], idx+offset, CI95[,2], lwd=1, col=col.pts)         
        points(idx+offset, coefs, pch=20, cex=cex.pts, col=col.pts)        
      }
      else{
        plot(c(idx+k,idx-k), c(CI95[,1],CI95[,2]), type="n",                                     
          axes=F, main=main, xlab=xlab, ylab=ylab,...) 
        if(plot){
          if (v.axis){
            axis(2)
          }
          if (h.axis){
            axis(1, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
              lty=0, cex.axis=cex.var)  
          }
          abline(h=0, lty=2)                                                 
          segments (idx, CI50[,1], idx, CI50[,2], lwd=2, col=col.pts)
          segments (idx, CI95[,1], idx, CI95[,2], lwd=1, col=col.pts)         
          points(idx, coefs, pch=20, cex=cex.pts, col=col.pts)
        }
      }
    }
  }
}
)
    


setMethod("coefplot", signature(object = "polr"), 
    function(object, varnames=NULL,...)
    {
    # collect informations
    coefs <- summary(object)$coef[,1]
    sds <- summary(object)$coef[,2]
    ifelse(is.null(varnames), varnames <- names(coefs), 
        varnames <- varnames)

    # plotting
    coefplot(coefs, sds, 
        varnames=varnames, ...)
    }
)  
