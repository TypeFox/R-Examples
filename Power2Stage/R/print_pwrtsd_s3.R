# ---------------------------------------------------------------
# S3 print method for class pwrtsd
# ---------------------------------------------------------------

print.pwrtsd <- function(x, ...)
{
  if(is.null(x$design)) x$design="2x2 crossover"
  cat("TSD with", x$design, "\n")
  
  # power.2stage.GS() has no method argument
  if(is.null(x$method)){
    cat(" non-adaptive group sequential with\n")
    cat("alpha (s1/s2) =", x$alpha[1], x$alpha[2], "\n")
    if(x$fCrange[1L] >0 & is.finite(x$fCrange[2L])){
      cat("Futility criterion ", x$fCrit," outside ", x$fCrange[1L], " ... ",
          x$fCrange[2L], "\n", sep="")
    } else {
      cat("No futility criterion\n")
    }
    cat("BE acceptance range = ", x$theta1," ... ", x$theta2,"\n\n", sep="")
    cat("CV= ", x$CV,"; n(s1, s2)= ", x$n[1]," " , x$n[2], "\n", sep="")
    
    # no sample size distribution in results here
    .print_results(x)

    return(invisible(x))
  }
  
  # power.2stage.ssr() has method="SSR" in return
  if(tolower(x$method)=="ssr"){
    if(x$blind) blinded <- "blinded " else blinded <- ""
    cat(" with ", blinded, "sample size re-estimation (only)\n", sep="")
    cat("Nominal alpha =", x$alpha, "\n")
    cat("Sample size est. based on power calculated via", 
        .pmethod_nice(x$pmethod), "\n")
    if(!x$usePE) {
      cat("with CV1, GMR =", x$GMR,"and targetpower =", x$targetpower,"\n")
    } else {
      cat("with CV1, PE1 and targetpower =", x$targetpower,"\n")
    }
    if(is.finite(x$max.n)){
      cat("Maximum sample size max.n = ", x$max.n,"\n", sep="")
    }
    if(x$min.n>0){
      cat("Minimum sample size min.n = ", x$min.n,"\n", sep="")
    }
    cat("BE acceptance range = ", x$theta1," ... ", x$theta2,"\n\n", sep="")
    cat("CV= ", x$CV,"; n(stage 1) = ", x$n1,"\n", sep="")
    
    .print_results(x)
    
    return(invisible(x))
    
  } 
  
  #all other functions
  cat("Method ", x$method,":", sep="")
  if (x$method=="C") cat(" alpha0 = ", x$alpha0, ",",sep="")
  cat(" alpha (s1/s2) =", x$alpha[1], x$alpha[2], "\n")
  
  # power.2stage.KM()
  if(!is.null(x$modified)){
    if(x$modified=="KM") {
      cat("Modification(s) according to Karalis/Macheras:\n")
      cat("- PE and mse of stage 1 in power steps AND sample size est. used\n") 
      cat("- Futility criterion PE outside 0.8 ... 1.25\n") 
    }
    if(x$modified=="DL") {
      cat("Modification according to DL:\n")
      cat("- PE and mse of stage 1 in sample size est. used if PE1 outside\n") 
      cat("  GMR ... 1/GMR, else GMR and mse1 is used\n") 
    }
  }

  # power.2stage.p()
  if(!is.null(x$test)){
    cat("CI's based on", ifelse(x$test=="welch","Welch's t-test", x$test),"\n")
  }
  if(is.null(x$powerstep)){
    cat("Target power in power monitoring and sample size est. = ", 
        x$targetpower,"\n",sep="")
    .print_pmethod(x$pmethod)
  } else {
    # only function power.2stage.fC() has argument powerstep
    if (x$powerstep) cat("Interim power monitoring step included\n") else 
      cat("No interim power monitoring step used\n")
    if (x$powerstep){
      cat("Target power in power monitoring and sample size est. = ", 
          x$targetpower,"\n",sep="")
    } else {
      cat("Target power in sample size est. = ", x$targetpower,"\n",sep="")
    }  
    .print_pmethod(x$pmethod)
  }

  # usePE is not in all results 
  if(!is.null(x$usePE)){
    if(x$usePE) cat("CV1 and PE1 in sample size est. used\n") else {
      cat("CV1 and GMR =", x$GMR, "in sample size est. used\n")}
  } 
  # futility criterion Nmax is not in all results
  if(!is.null(x$Nmax)){
    if(is.finite(x$Nmax)) {
      cat("Futility criterion Nmax = ",x$Nmax,"\n", sep="")
    } else {
      cat("No futility criterion\n")
    }
  }
  # max.n is not always present
  if(!is.null(x$max.n)){
    if(is.finite(x$max.n)){
      cat("Maximum sample size max.n = ", x$max.n,"\n", sep="")
    }
  }
  
  # futility criterion w.r.t. PE or CI
  if(!is.null(x$fCrit)){
    if(x$fCrange[1L] >0 & is.finite(x$fCrange[2L])){
      fCrit <- x$fCrit
      if(fCrit=="CI") fCrit <- paste0(100*(1-2*x$alpha0), "% CI")
      cat("Futility criterion ", fCrit," outside ", x$fCrange[1L], " ... ",
          x$fCrange[2L], "\n", sep="")
    } else {
      cat("No futility criterion\n")
    }
  }
  if(!is.null(x$min.n2)){
    if(x$min.n2>0) cat("Minimum sample size in stage 2 =", x$min.n2, "\n")
  }
  # --- BE acceptance range
  cat("BE acceptance range = ", x$theta1," ... ", x$theta2,"\n\n", sep="")
  
  # ---CV, n1 and GMR (GMR only if usePE=F)
  if (length(x$CV)==2){
    cat("CVs = "); cat(x$CV, sep=", ")
  } else {
    cat("CV = ", x$CV, sep="")
  }
  if (length(x$n1)==1){
    cat("; n(stage 1) = ", x$n1, sep="")
  } else {
    #this is the case for power.2stage.p only
    cat("; ntot(stage 1) = ", sum(x$n1), " (nT, nR = ", x$n1[1], ", ", x$n1[2],")", sep="")
  }
  
  # power.2stage.KM() has no GMR
  # if usePE=T then also no output of GMR
  if(!is.null(x$usePE)) {
    if(x$usePE==TRUE) x$GMR <- NULL
  } 
  if(!is.null(x$GMR)) cat("; GMR=", x$GMR, "\n") else cat("\n")
  
  .print_results(x)

  return(invisible(x))
  
} # end function

# nice power calc. method
.print_pmethod <- function(pmethod){
  cat("Power calculation via", .pmethod_nice(pmethod),"\n")
}

.pmethod_nice <- function(pmethod)
{
  if(pmethod=="exact")   ret <- "exact method"
  if(pmethod=="nct")     ret <- "non-central t approx."
  if(pmethod=="shifted") ret <- "shifted central t approx."
  if(pmethod=="ls")      ret <- "large sample normal approx."
  if(pmethod=="exp")     ret <- "'expected power'"
  ret
}

.print_results <- function(x)
{
  # --- results part
  cat("\n", x$nsims," sims at theta0 = ", x$theta0, sep="")
  if(x$theta0<=x$theta1 | x$theta0>=x$theta2) cat(" (p(BE)='alpha').\n") else { 
    cat(" (p(BE)='power').\n")}
  cat("p(BE)    = ", x$pBE,"\n", sep="")
  # power.2stage.ssr() has no pBE_s1
  if(!is.null(x$pBE_s1)){
    cat("p(BE) s1 = ", x$pBE_s1,"\n", sep="")
  }
  cat("Studies in stage 2 = ", round(x$pct_s2,2), "%\n", sep="")
  # power.2stage.GS() has no sample size distribution at all
  # we assume if nmean is given we have also the range
  if(!is.null(x$nmean)){
    cat("\nDistribution of n(total)\n")
    cat("- mean (range) = ", round(x$nmean,1)," (", x$nrange[1]," ... ",
        x$nrange[2],")\n", sep="")
    # percentiles are not in all results
    if(!is.null(x$nperc)){
      cat("- percentiles\n")
      print(x$nperc)
    }
  }
  cat("\n")
}
