# -------------------------------------------------------------------
# power analysis for NTID according to FDA method
# Author: D. Labes after pa.ABE() from Helmut Schuetz 
# Only design="2x2x4" according to the Warfarin guidance
# -------------------------------------------------------------------
pa.NTIDFDA <- function(CV, theta0=0.975, targetpower=0.8, minpower=0.7,  ...) 
{ 
  # Rversion must be >=3.1.0 for the uniroot call with argument extendInt
  Rver <- paste0(R.Version()$major, ".", R.Version()$minor)
  
  # functions to use with uniroot
  pwrCV  <- function(x, ...) {
    power.NTIDFDA(CV=x, ...) - minpower
  } 
  pwrGMR <- function(x, ...) {
    power.NTIDFDA(theta0=x, ...) - minpower
  } 
  # to avoid reprogramming of Helmuts code
  GMR <- theta0
  
  if (length(CV)>1) stop("Only CVwT=CVwR=CV implemented. CV has to be a scalar.")
  if (CV<0) {
    CV <- -CV
    message("Negative CV changed to ",CV,".")
  }
  if (targetpower>=1 | minpower>=1 | targetpower<=0 | minpower<=0 ) 
    stop("Power values have to be within 0...1") 
  if(minpower>=targetpower) stop("Minimum acceptable power must < than target.")
  if (Rver<"3.1.0"){
    if (minpower < 0.5) stop("Minimum acceptable power must be >=0.5.")
    if (targetpower < 0.5) stop("Target power must be >=0.5.")
  } else {
    if (minpower <= 0.5) 
      message("Note: Minimum acceptable power <=0.5 doesn't make much sense.")
    if (targetpower <= 0.5) 
      message("Note: Target power <=0.5 doesn't make much sense.")
  }
  
  #if(minpower < 0.05) stop("Minimum acceptable power must not be <0.05!")
  # otherwise uniroot() throws an error says Helmut
  
  # power.NTIDFDA has only full replicate design implemented
  design <- "2x2x4"
  dno     <- .design.no(design)
  if (is.na(dno)) stop("Design,", design, " not implemented")
  d.props <- .design.props(dno)
  seqs    <- d.props$steps
  
  res  <- sampleN.NTIDFDA(CV=CV, theta0=GMR, targetpower=targetpower, 
                          print=FALSE, details=FALSE, ...)
  n.est   <- res[1, "Sample size"   ]
  pwr.est <- res[1, "Achieved power"]
  # don't allow below 12 subjects
  # may not necessary, have upto now 14 seen as min.
  if(n.est<12){
    n.est <- 12
    pwr.est <- power.NTIDFDA(CV=CV, n=n.est, ...)
    res[,"Sample size"] <- n.est
    res[,"Achieved power"] <- pwr.est
  }
  ########################################
  # max. CV for minimum acceptable power #
  ########################################
  # TODO: for small CV's and theta0^=1 the search interval may be not sufficient
  # to get a solution in uniroot
  upr <- max(CV*10, 1/CV)
  CV.max <- NA
  # for the try code see Hadley Wickhams "Advanced R"
  try(CV.max <- uniroot(pwrCV,  c(CV, upr), tol=1e-7, n=n.est, 
                        theta0=GMR, ...)$root, silent = TRUE)
  CV.min <- NA
  try(CV.min <- uniroot(pwrCV,  c(0.0001, CV), tol=1e-7, n=n.est, 
                        theta0=GMR, ...)$root, silent = TRUE)
  if (is.na(CV.min)) CV.min <- CV*2/3 else CV.min <- as.numeric(CV.min)
  if (is.na(CV.max)) CV.max <- 0.5    else CV.max <- as.numeric(CV.max)
  # points for plotting
  seg    <- 60; 
  # DL: simplified code of Helmut
  # do we need here different stepsizes
  CVs    <- c(seq(CV.min, CV, length.out=seg*.25), 
              seq(CV, CV.max, length.out=seg*.75)[-1])
  
  # dimension vector properly in advance
  pBECV  <- vector("numeric", length=length(CVs))
  # replace next eventually with apply/sapply
  for(j in seq_along(CVs)) {
    pBECV[j] <- power.NTIDFDA(CV=CVs[j], n=n.est, theta0=GMR, ...)
  }
  
  ######################################
  # min. GMR for minimum accept. power #
  ######################################
  seg      <- 50
  if(GMR <= 1) interval <- c(0.8, GMR) else interval <- c(GMR, 1.25)
  GMR.min  <- uniroot(pwrGMR, interval, tol=1e-7, n=n.est, CV=CV, ...)$root
  GMRs     <- seq(GMR.min, GMR, length.out=seg)
  pBEGMR   <- vector("numeric", length=length(GMRs))
  # replace next loop eventually with apply/sapply/vapply
  for(j in seq_along(GMRs)) {
    pBEGMR[j] <- power.NTIDFDA(CV=CV, n=n.est, theta0=GMRs[j], ...)
  }
  
  ####################################
  # min. n for minimum accept. power #
  # workaround, since uniroot() does # what does he mean here?
  # not accept two vectors as limits #
  ####################################
  nstep <- 1
  if (n.est>200) nstep <- 6
  if (n.est>500) nstep <- 10
  if (n.est>1000) nstep <- 50
  if (n.est>5000) nstep <- 100
  
  Ns   <- seq(n.est, 12, by=-nstep) # don't drop below 12 subjects
  if (n.est==12) Ns <- seq(n.est, 8)
  nNs  <-length(Ns)
  j    <- 0
  pwrN <- pwr.est
  # may it be that j grows greater than length(Ns)?
  n.min <- NULL; pBEn <- NULL
  # vector of n in sequence groups
  n <- vector("numeric", length=seqs)
  ni <- 1:seqs
  while(pwrN >= minpower & j<nNs){
    j <- j+1
    # next is done in power.NTIDFDA(), but left here
    n[-seqs] <- diff(floor(Ns[j]*ni/seqs))
    n[seqs]  <- Ns[j] -sum(n[-seqs])
    # debug prints
    #cat("Ns[j]",Ns[j],"\n")
    #print(n)
    pwrN <- power.NTIDFDA(CV=CV, n=n, theta0=GMR, ...)
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
             method="RSABE NTID",
             regulator="FDA"
             )
  class(ret) <- c("pwrA", class(ret))

  return(ret)

}