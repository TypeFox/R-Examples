# -----------------------------------------------------------
# power analysis 
# Author: Helmut Schuetz 
# with some modifications by D. Labes to adapt to PowerTOST
# infrastructure and namings
# -----------------------------------------------------------
pa.ABE <- function(CV, theta0=0.95, targetpower=0.8, minpower=0.7, 
                   design="2x2", ...) 
{ # Rversion must be >=3.1.0 for the uniroot call with argument extendInt
  Rver <- paste0(R.Version()$major, ".", R.Version()$minor)
  
  # functions to use with uniroot
  pwrCV  <- function(x, ...) {
    power.TOST(CV=x, ...) - minpower
  } 
  pwrGMR <- function(x, ...) {
    power.TOST(theta0=x, ...) - minpower
  } 
  # to avoid reprogramming of Helmuts code
  GMR <- theta0
  
  if (targetpower>=1 | minpower>=1) stop("Power values have to be within 0...1") 
  if (targetpower<=0 | minpower<=0) stop("Power values have to be within 0...1") 
  if (minpower>=targetpower) stop("Minimum acceptable power must < than target")
  if (Rver<"3.1.0"){
    if (minpower < 0.5) stop("Minimum acceptable power must be >=0.5.")
    if (targetpower < 0.5) stop("Target power must be >=0.5.")
    
  } else {
    if (minpower <= 0.5) 
      message("Note: Minimum acceptable power <=0.5 doesn't make much sense.")
    if (targetpower <= 0.5) 
      message("Note: Target power <=0.5 doesn't make much sense.")
  }
  
  if (CV<0) {
    CV <- -CV
    message("Negative CV changed to ",CV,".")
  }
  
  dno     <- .design.no(design)
  if (is.na(dno)) stop("Design,", design, " not implemented")
  d.props <- .design.props(dno)
  seqs    <- d.props$steps
  
  res  <- sampleN.TOST(CV=CV, theta0=GMR, targetpower=targetpower, 
                       design=design, print=FALSE, details=FALSE, ...)
  n.est   <- res[1, "Sample size"   ]
  pwr.est <- res[1, "Achieved power"]
  
  # don't allow below 12 subjects
  if(n.est<12){
    n.est <- 12
    pwr.est <- power.TOST(CV=CV, n=n.est, design=design, ...)
    res[,"Sample size"] <- n.est
    res[,"Achieved power"] <- pwr.est
  }
  # points for plotting
  seg     <- 50; s <- seq(seg-1)
  
  ########################################
  # max. CV for minimum acceptable power #
  ########################################
  if(Rver<"3.1.0"){
    CV.max <- uniroot(pwrCV,  c(CV, 10*CV), tol=1e-7, 
                      n=n.est, theta0=GMR, design=design, ...)$root
  } else {
    # argument extentInt is available from 3.1.0 on
    CV.max <- uniroot(pwrCV,  c(CV, 10*CV), tol=1e-7, extendInt ="downX",
                      n=n.est, theta0=GMR, design=design, ...)$root
  }
  CVs    <- seq(CV, CV.max, length.out=seg)
  # dimension properly in advance
  pBECV  <- vector("numeric", length=length(CVs))
  # 1:length(CVs) is R-Inferno, also growing objects in steps
  # replace eventually with apply/sapply
  for(j in seq_along(CVs)) {
    pBECV[j] <- power.TOST(CV=CVs[j], n=n.est, theta0=GMR, design=design, ...)
  }
  
  ######################################
  # min. GMR for minimum accept. power #
  ######################################
                          # original was here (0.8, 1)
  if(GMR <= 1) interval <- c(0.8, GMR) else interval <- c(GMR, 1.25)
  # extending interval for extrem cases
  if(GMR <= 1) updown   <- "upX" else updown <- "downX"
  if(Rver<"3.1.0"){
    GMR.min  <- uniroot(pwrGMR, interval, tol=1e-7, 
                        n=n.est, CV=CV, design=design, ...)$root
  } else {
    # extentInt is available from 3.1.0 on
    GMR.min  <- uniroot(pwrGMR, interval, tol=1e-7, extendInt=updown,
                        n=n.est, CV=CV, design=design, ...)$root
  }  
  GMRs     <- seq(GMR.min, GMR, length.out=seg)
  pBEGMR   <- vector("numeric", length=length(GMRs))
  # replace eventually with apply/sapply
  for(j in seq_along(GMRs)) {
    pBEGMR[j] <- power.TOST(CV=CV, n=n.est, theta0=GMRs[j], design=design, ...)
  }
  
  ####################################
  # min. n for minimum accept. power #
  # workaround, since uniroot() does #
  # not accept two vectors as limits #
  ####################################
  #Ns <- seq(n.est, 12, by=-1) # don't drop below 12 subjects
  Ns  <- seq(n.est, 12)
  if(n.est==12) Ns <- seq(n.est, 2*seqs)
  nNs <-length(Ns)
  j   <- 0
  pwrN <- pwr.est
  # may it be that j grows greater than length(Ns)?
  n.min <- NULL; pBEn <- NULL
  n <- vector("numeric", length=seqs)
  ni <- 1:seqs
  while(pwrN >= minpower & j<nNs){
    j <- j+1
    n[-seqs] <- diff(floor(Ns[j]*ni/seqs))
    n[seqs]  <- Ns[j] -sum(n[-seqs])
    pwrN <- power.TOST(CV=CV, n=n, theta0=GMR, design=design, ...)
    if(pwrN >= minpower) {
      n.min <- c(n.min, sum(n))
      pBEn <- c(pBEn, pwrN)
    } else {
      break
    }
  }
  
  # plots are now contained in the S3 method plot

  # return what?
  ret <-list(plan=res, 
             paCV=data.frame(CV=CVs, pwr=pBECV),
             paGMR=data.frame(theta0=GMRs, pwr=pBEGMR),
             paN=data.frame(N=n.min, pwr=pBEn),
             minpower=minpower,
             method="ABE"
             )
  class(ret) <- c("pwrA", class(ret))

  return(ret)

}