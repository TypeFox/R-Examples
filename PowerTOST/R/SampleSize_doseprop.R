#source(helper_dp.R)
# -----------------------------------------------------------------------
# large sample estimate of sample size
# -----------------------------------------------------------------------
sampleN0.dp <- function(alpha=0.05, CV, doses, targetpower=0.8, beta0=1, 
                        theta1=0.8, theta2=1/theta1, design="crossover", 
                        dm=NULL, CVb)
{
  s2   <- CV2mse(CV)
  
  rd   <- max(doses)/min(doses)
  
  grps <- length(doses)
  prds <- ifelse(design=="parallel", 1, grps)
  prds <- ifelse(design=="IBD",ncol(dm), prds)
  grps <- ifelse(design=="IBD",nrow(dm), grps)
  # all ni=1
  n <- rep.int(1, times=grps)
  # corrected sum of squares of log doses of design
  css  <- .css3(doses, design, dm, n, s02=CV2mse(CV), omega2=CV2mse(CVb))
  # lower / upper acceptance range of slope
  bl   <- 1+log(theta1)/log(rd)
  bu   <- 1+log(theta2)/log(rd)
  
  if (beta0<=bl | beta0>=bu){
    warning("beta0 must lie between ",signif(bl,5)," and ",
            signif(bu,5),". NaN returned.")
    return(NaN)
  } 
  
  # quantil for type2 error: U(power) or U(1-beta/2)
  if (beta0!=1) Ub <- qnorm(targetpower) else Ub <- qnorm(1-(1-targetpower)/2)
  # see Julious, Tan & Machin formula for beta0=1
  # see also Sethuraman et al.
  n <- s2*(Ub+qnorm(1-alpha))^2/css
  if (beta0>=1) n <- n/(beta0-bu)^2 else n<- n/(beta0-bl)^2
  n <- grps*round(n,0)
  return(n)
}

# -----------------------------------------------------------------------
# sample size for dose prop., power model
# -----------------------------------------------------------------------
sampleN.dp <- function(alpha=0.05, CV, doses, targetpower=0.8, beta0, 
                       theta1=0.8, theta2=1/theta1, 
                       design=c("crossover", "parallel", "IBD"), dm=NULL,
                       CVb, print=TRUE, details=FALSE, imax=100)
{
  design <- match.arg(design)
  
  grps <- length(doses) # dose groups and periods in case of crossover
  if (grps<=1) stop("At least two doses have to be given.")
  if (design=="IBD"){
    # check
    if(!is.matrix(dm)) stop("Design matrix must be given.")
    grps <- nrow(dm) # sequence groups
    if(missing(CVb)) CVb <- 2*CV
  } else {
    if(missing(CVb)) CVb <- 0
  }
  
  if (CV<=0) stop("CV must be greater then zero.")
  s2   <- CV2mse(CV)
  
  # acceptance range for slope
  rd   <- max(doses)/min(doses)
  bl   <- 1+log(theta1)/log(rd)
  bu   <- 1+log(theta2)/log(rd)
  
  if (missing(beta0)) beta0 <- 1+log(0.95)/log(rd)
  if (beta0<=0) stop("beta0 must be greater then zero.")
  if (beta0<=bl | beta0>=bu) stop("beta0 must lie between ",signif(bl,5)," and ",
                                   signif(bu,5),".")
   
  if (theta1<=0 | theta2<=0) stop("theta1/theta2 must be greater then zero.")
  
  if (print){
    cat("\n++++ Dose proportionality study, power model ++++\n")
    cat("            Sample size estimation\n")
    cat("-------------------------------------------------\n")
    cat("Study design: ",design, sep="")
    descomm <- paste0(" (",grps,"x",grps," Latin square)")
    if (design=="parallel"){
      descomm <- paste0(" (",grps," groups)")
    }
    if (design=="IBD"){
      descomm <- paste0(" (",length(doses),"x",grps,"x",ncol(dm),")")
    }
    cat(descomm,"\n")
    cat("alpha = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("Equivalence margins of R(dnm) =",theta1,"...", theta2,"\n")
    cat("Doses = "); cat(doses,"\n")
    cat("Null (true) slope = ",beta0,", CV = ", CV, sep="")
    if(design=="IBD"){
      cat(", CVb = ", CVb, "\n", sep="")
    } else {
      cat("\n")
    }
    cat("Slope acceptance range =",signif(bl,5),"...", signif(bu,5),"\n")
  }
  # start sample size search
  n <- sampleN0.dp(alpha=alpha, CV=CV, doses=doses, targetpower=targetpower, 
                   beta0=beta0, theta1=theta1, theta2=theta2, design=design,
                   dm=dm, CVb=CVb)
  # here for small CV < 0.1 n=0 may result
  # is this reasonable?
  if (n<grps)  n <- grps
  if (grps==2 & n<4) n <- 4  # else df=0 may result for "crossover" or "parallel"
  pwr <- power.dp(alpha=alpha, CV=CV, doses=doses, n=n, beta0=beta0, 
                  theta1=theta1, theta2=theta2, design=design, dm=dm, CVb=CVb)
  if (details) {
    cat("\nSample size search (ntotal)\n")
    cat(" n     power\n")
    # do not print first too high
    # this is for cases with only one step-down and than step up
    if (pwr<=targetpower) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
  }
  iter <- 0
  nmin <- 3
  while(pwr>targetpower){
    if (n<=nmin) { 
      if (details & iter==0) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
      break
    }
    # step down
    n   <- n-grps
    pwr <- power.dp(alpha=alpha, CV=CV, doses=doses, n=n, beta0=beta0, 
                    theta1=theta1, theta2=theta2, design=design, dm=dm, CVb=CVb)
    iter <- iter+1
    if (details) cat( n," ", formatC(pwr, digits=6),"\n")
    if (iter>imax) break  
  }
  while(pwr<targetpower){
    # step up
    n   <- n+grps
    pwr <- power.dp(alpha=alpha, CV=CV, doses=doses, n=n, beta0=beta0, 
                    theta1=theta1, theta2=theta2, design=design, dm=dm, CVb=CVb)
    if (details) cat( n," ", formatC(pwr, digits=6),"\n")
    iter <- iter+1
    if (iter>imax) break  
  }
  if(print && !details){
    cat("\nSample size (total)\n")
    cat(" n     power\n")
    cat( n," ", formatC(pwr, digits=6, format="f"), "\n")
  }
  #return results as data.frame
  res <- data.frame(design=design, alpha=alpha, CV=CV, CVb=CVb,
                    doses=paste(doses, collapse=", "), beta0=beta0, 
                    theta1=theta1, theta2=theta2, n=n, power=pwr, 
                    targetpower=targetpower)
  names(res) <-c("Design", "alpha", "CV", "CVb", "doses", "beta0", "theta1", "theta2",
                 "Sample size", "Achieved power", "Target power")
  if (design!="IBD") res$CVb <- NULL
  
  if (print) return(invisible(res)) else return(res)
  
} #end function
  