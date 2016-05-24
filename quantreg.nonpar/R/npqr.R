npqr <-
function(formula, data, basis=NULL, var, taus=c(.25,.5,.75), 
                           print.taus=NULL, B=200, nderivs=1, average=T, 
                           load=NULL, alpha=0.05, process="pivotal", rearrange=F, rearrange.vars="quantile",
                           uniform=F, se="unconditional", printOutput=T, method="fn"){

  if (is.null(print.taus)) { # If taus to be printed are not supplied, set to all taus
    print.taus <- taus
  }
  if (any(!(print.taus %in% taus))) {
    stop("print.taus must be a subset of taus")
  }

  # Generate main output
  if (process=="pivotal") {
    V    <- pivotal(data=data, B=B, taus=taus, formula=formula, basis=basis, alpha=alpha, var=var, load=load,
                rearrange=rearrange, rearrange.vars=rearrange.vars, uniform=uniform, se=se, average=average, nderivs=nderivs, method=method)
  } else if (process=="gaussian") {
    V    <- gaus(data=data, B=B, taus=taus, formula=formula, basis=basis, alpha=alpha, var=var, load=load,
                rearrange=rearrange, rearrange.vars=rearrange.vars, uniform=uniform, se=se, average=average, nderivs=nderivs, method=method)
  } else if (process=="wbootstrap"){
    V    <- wbootstrap(data=data, B=B, taus=taus, formula=formula, basis=basis, alpha=alpha, var=var, load=load, 
                       rearrange=rearrange, rearrange.vars=rearrange.vars, uniform=uniform, average=average, nderivs=nderivs, method=method) 
  } else if (process=="gbootstrap"){
    V    <- gbootstrap(data=data, B=B, taus=taus, formula=formula, basis=basis, alpha=alpha, var=var, load=load, 
                       rearrange=rearrange, rearrange.vars=rearrange.vars, uniform=uniform, average=average, nderivs=nderivs, method=method)
  } else if (process=="none"){
    V    <- no.process(data=data, taus=taus, formula=formula, basis=basis, var=var, load=load, 
                       rearrange=rearrange, rearrange.vars=rearrange.vars, average=average, nderivs=nderivs, method=method)
  }

  length.taus <- length(taus)
  loadLength <- dim(V$load)[1]
  
  output <- list(CI=NULL, CI.oneSided=NULL, point.est=NULL, std.error=NULL, pvalues=NULL, taus=NULL, coefficients=NULL, var.unique=NULL, load=NULL)
  if (process!="none"){
    output$CI <- V$CI
    output$CI.oneSided <- V$CI.oneSided
    output$std.error <- V$std.error
  }
  output$taus <- taus
  output$coefficients <- vector(mode="list", length=length.taus)
  for (i in 1:length.taus) {
    output$coefficients[[i]] <- V$qfits[[i]]$coef
  }
  
  output$var.unique <- V$var.unique
  output$point.est <- V$point.est
  output$pvalues <- V$pvalues
  output$load <- V$load

  # Print output

  nObs <- length(V$qfits[[1]]$residuals)

  if(printOutput==T){
    cat("\n\n")
    if (process=="pivotal"){
      cat(format("Method: pivotal"),"\n")
    } else if (process=="gaussian"){
      cat(format("Method: gaussian"),"\n")
    } else if (process=="wbootstrap"){
      cat(format("Method: weighted bootstrap"),"\n")
    } else if (process=="gbootstrap"){
      cat(format("Method: gradient bootstrap"),"\n")
    } else if (process=="none"){
      cat(format("Method: no inference"), "\n")
    }
    if (average==T & (process=="gaussian" | process=="pivotal") & se=="unconditional"){
      cat(format("Standard error estimation is unconditional."),"\n\n")
    } else if (average==T & (process=="gaussian" | process=="pivotal") & se!="unconditional"){
      cat(format("Standard error estimation is conditional."),"\n\n")
    } else {
      cat("\n")
    }
    cat(format("No. of obs.:",width=50),nObs,"\n")
    if (process!="none"){
      cat(format("No. of simulations or bootstrap repetitions:",width=50),B,"\n\n\n")
      cat(paste(rep("-",80),collapse=""),"\n")
      cat(format("Quantile Estimates & Inference", width=80,justify = "centre"),"\n")
      cat(paste(rep("-",80),collapse=""),"\n")
      if (uniform==T){
        cat(format("",width=10),format("Point",width=10,justify="centre"), format("Standard",width=10,justify="centre"), format("Uniform",width=20,justify="centre"), format(paste("One-sided ",(1-alpha)*100, "% CIs",sep=""),width=20,justify="centre"),"\n")
      } else {
        cat(format("",width=10),format("Point",width=10,justify="centre"), format("Standard",width=10,justify="centre"), format("Pointwise",width=20,justify="centre"), format(paste("One-sided ",(1-alpha)*100, "% CIs",sep=""),width=20,justify="centre"),"\n")
      }
      cat(format("Quantile",width=10,justify="right"),format("Estimate",width=10,justify = "centre"),format("Error",width=10,justify="centre"), format(paste((1-alpha)*100,"% Conf.Interval",sep=""),width=20,justify="centre"),format("Lower",width=10,justify = "centre"),format("Upper",width=10,justify = "centre"),"\n")    	
      for (j in 1:loadLength){
        for(i in 1:length(taus)){
          if (taus[i] %in% print.taus){
            cat(format(taus[i],width=10,justify="left"),format(V$load[j,] %*% V$qfits[[i]]$coef,justify="centre",width=10,digits=4),format(output$std.error[j,i],justify="centre",width=10,digits=4), format(output$CI[j,i,1],justify="centre",width=10,digits=4),format(output$CI[j,i,2],justify="centre",width=10,digits=4),format(output$CI.oneSided[j,i,1],justify="centre",width=10,digits=4),format(output$CI.oneSided[j,i,2],justify="centre",width=10,digits=4),"\n",sep="")
          }
        }
        cat("\n")
        cat(paste(rep("-",80),collapse=""),"\n")
      }
      cat(format("Hypothesis Testing", width=80,justify = "centre"),"\n")
      cat(paste(rep("-",80),collapse=""),"\n")
      cat(format("Null Hypothesis:", width=40, justify = "left"), format("p-value", width=40, justify="right"), "\n")
      cat(paste(rep("=",80),collapse=""), "\n")
      cat(format("theta(tau, w) <= 0 for all tau, w", width=40, justify= "left"), format(output$pvalues[1], width=40, justify="right", digits=4), "  \n")
      cat(format("theta(tau, w) >= 0 for all tau, w", width=40, justify= "left"), format(output$pvalues[2], width=40, justify="right", digits=4), "  \n")
      cat(format("theta(tau, w) = 0 for all tau, w", width=40, justify= "left"), format(output$pvalues[3], width=40, justify="right", digits=4), "  \n")
      cat(paste(rep("=",80),collapse=""))
    } else {
      cat(paste(rep("-",80),collapse=""),"\n")
      cat(format("Quantile Estimates", width=80,justify = "centre"),"\n")
      cat(paste(rep("-",80),collapse=""),"\n")
      cat(format("Quantile",width=30,justify="right"),format("Point Estimate",width=20,justify="right"),"\n")
      for (j in 1:loadLength){
        for(i in 1:length(taus)){
          if (taus[i] %in% print.taus){
            cat(format(taus[i],width=30,justify="right"),format(V$load[j,] %*% V$qfits[[i]]$coef,justify="right",width=20,digits=4),"\n",sep="")
          }
        }
        cat("\n")
        cat(paste(rep("-",80),collapse=""),"\n")
      }
    }
  }

 return(output)
}
