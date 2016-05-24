#require(PowerTOST)
#-------------------------------------------------------------------------
# power analysis of a sample size plan for scaled ABE (EMA & FDA)
# Input area at the end of the code.
# coded originally by Helmut Schuetz
# adapted to PowerTOST infra structure by D. Labes
#-------------------------------------------------------------------------
pa.scABE   <- function(CV, theta0=0.9, targetpower=0.8, minpower=0.7, 
                       design=c("2x3x3", "2x2x4", "2x2x3"), 
                       regulator=c("EMA", "ANVISA", "FDA"), ...) 
{ # Rversion must be >=3.1.0 for the uniroot call with argument extendInt
  Rver <- paste0(R.Version()$major, ".", R.Version()$minor)
  
  if (length(CV)>1) stop("Only CVwT=CVwR=CV implemented. CV has to be a scalar.")
  if (CV<0) {
    CV <- -CV
    message("Negative CV changed to ",CV,".")
  }
  
  if (targetpower>=1 | minpower>=1 | targetpower<=0 | minpower<=0 ) 
    stop("Power values have to be within 0...1") 
  if(minpower>=targetpower) stop("Minimum acceptable power must be < than target")
  
  if (Rver<"3.1.0"){
    if (minpower < 0.5) stop("Minimum acceptable power must be >=0.5.")
    if (targetpower < 0.5) stop("Target power must be >=0.5.")
  } else {
    if (minpower <= 0.5) 
      message("Note: Minimum acceptable power <=0.5 doesn't make much sense.")
    if (targetpower <= 0.5) 
      message("Note: Target power <=0.5 doesn't make much sense.")
  }
  
  # check regulator, check design
  regulator <- toupper(regulator)
  reg       <- match.arg(regulator)
  design    <- match.arg(design)
  
  # to avoid recoding of Helmut's code
  GMR <- theta0
  
  # functions to be used with uniroot
  # reg is visible as long as these functions are def. within pa.scABE()
	pwrCV <- function(x, ...) 
  {
    if (reg=="FDA") {
      power.RSABE(CV=x, ...) - minpower
    } else {
      power.scABEL(CV=x, ...)- minpower
    }
	}
  pwrGMR <- function(x, ...) 
  {
    if (reg == "FDA"){
      power.RSABE(theta0=x, ...) - minpower
    } else {
      power.scABEL(theta0=x, ...) - minpower
    }
  }
  
  if(reg == "FDA") {
    res <- sampleN.RSABE(CV=CV, theta0=GMR, targetpower=targetpower, 
                         design=design, print=FALSE, details=FALSE, ...)
	} else {
	  res <- sampleN.scABEL(CV=CV, theta0=GMR, targetpower=targetpower, 
	                        design=design, regulator=reg,
                          print=FALSE, details=FALSE, ...)
	}
	n.est   <- res[1, "Sample size"   ]
	pwr.est <- res[1, "Achieved power"]
  
  # don't allow below 12 subjects
  if(n.est<12){
    n.est <- 12
    if(reg == "FDA") {
      pwr.est <- power.RSABE(CV=CV, n=n.est, theta0=GMR,  design=design, ...)
    } else {
      pwr.est <- power.scABEL(CV=CV, n=n.est, theta0=GMR, design=design, 
                              regulator=reg, ...)
    }
    res[,"Sample size"] <- n.est
    res[,"Achieved power"] <- pwr.est
  }
  
  ########################################
	# max. CV for minimum acceptable power #
	########################################
  # TODO: test many cases. the form of the curves suggest that uniroot may fail!
  if (Rver<"3.1.0"){
    CV.max <- uniroot(pwrCV, c(CV, 30*CV), tol=1e-7,  
                      n=n.est, design=design, theta0=GMR, regulator=reg, ...)$root
    
  } else {
    CV.max <- uniroot(pwrCV, c(CV, 30*CV), tol=1e-7, extendInt ="downX", 
                      n=n.est, design=design, theta0=GMR, regulator=reg, ...)$root
  }
  # points for plotting
  seg    <- 75; 
  # what did he do here? D.L.
  CVs    <- c(seq(CV*2/3, CV, length.out=seg*.25), seq(CV, CV.max, length.out=seg*.75))
  CVs    <- c(head(CVs, floor(seg*0.25)), tail(CVs, seg*0.75)) # amateur! [Helmuts comment for himself :-))]
	pBECV  <- vector("numeric", length=length(CVs))
	for(j in seq_along(CVs)) {
		if(reg == "FDA") {
		  pBECV[j] <- power.RSABE(CV=CVs[j], n=n.est, theta0=GMR,  design=design, ...)
		} else {
		  pBECV[j] <- power.scABEL(CV=CVs[j], n=n.est, theta0=GMR, design=design,
                               regulator=reg, ...)
		}
	}
	
  ##################################################################
	# min. GMR for minimum accept. power, may also be max if GMR>1!  #
	##################################################################
	if(reg == "EMA" | reg=="ANVISA") { 
    # scABEL (rounded switch acc. to BE-GL and Q&A document)
		ifelse(CV <= 0.5, UL <- exp(0.76*CV2se(CV)), UL <- exp(0.76*CV2se(0.5)))
		if(CV <= 0.3) UL <- 1.25
    if (reg=="ANVISA"){
      if(CV <= 0.4) UL <- 1.25
    }
	} else {       
    # RSABE (exact switching cond.: 0.8925742...)
		ifelse(CV > 0.3, UL <- exp(log(1.25)/0.25*CV2se(CV)), UL <- 1.25)
	}
	ifelse(GMR <= 1, interval <- c(1/UL, 1), interval <- c(GMR, UL)) # scaled borders!
  # extending interval for extrem cases
  updown <- ifelse(GMR <= 1, "upX", "downX")
  
  if (Rver<"3.1.0"){
    GMR.min <- uniroot(pwrGMR, interval, tol=1e-7, 
                       n=n.est, CV=CV, design=design, regulator=reg, ...)$root
  } else {
    GMR.min <- uniroot(pwrGMR, interval, tol=1e-7, extendInt=updown,
                       n=n.est, CV=CV, design=design, regulator=reg, ...)$root
  }  
  seg     <- 50 # save some speed compared to Helmuts always 75
	GMRs    <- seq(GMR.min, GMR, length.out=seg)
	pBEGMR  <- vector("numeric", length=length(GMRs))
	for(j in seq_along(GMRs)) {
		if(reg == "FDA") {
		  pBEGMR[j] <- power.RSABE(CV=CV, n=n.est, theta0=GMRs[j], 
		                           design=design, ...)
		} else {
		  pBEGMR[j] <- power.scABEL(CV=CV, n=n.est, theta0=GMRs[j], 
		                            design=design, regulator=reg, ...)
		}
	}
	####################################
	# min. n for minimum accept. power #
	# workaround, since uniroot() does #
	# not accept two vectors as limits #
	####################################
	#Ns <- seq(n.est, 12, by=-1) # don't drop below 12 subjects (original)
  Ns    <- seq(n.est, 12)
  if(n.est==12) Ns <- seq(n.est, 6)
	pwrN  <- pwr.est
	n.min <- NULL; pBEn <- NULL
  #######################################################
	# THX to Detlew Labes for this part of the code!     #
	# See: http://forum.bebac.at/forum_entry.php?id=13375 #
	#######################################################
	# get # of sequences 
  # dont need this here since unbalancedness is handled
  # by power.scABEL() or power.RSABE() itself
# 	seqs <- known.designs()[known.designs()$design==design,"steps"]
# 	n    <- vector("numeric", length=seqs)
# 	ni   <- 1:seqs
	# may it be that j grows greater than length(Ns)? paranoia
	nNs <- length(Ns)
  j   <- 0
  while(pwrN >= minpower & j<nNs){
		j <- j+1
		# distribute total Ns to the sequence groups
		# n[-seqs] is n[1:(seqs-1)]
    # don't need this here, is handled in the power functions
# 		n[-seqs] <- diff(floor(Ns[j]*ni/seqs))
# 		n[seqs]  <- Ns[j] -sum(n[-seqs])
    # suppress messages regarding unbalanced designs
    suppressMessages(
  		if(reg == "FDA") {
  		  pwrN <- power.RSABE(CV=CV, n=Ns[j], theta0=GMR, design=design, ...)
  		} else {
  		  pwrN <- power.scABEL(CV=CV, n=Ns[j], theta0=GMR, design=design, 
                             regulator=reg, ...)
  		}
    )  
		if(pwrN >= minpower) {
			n.min <- c(n.min, Ns[j])
			pBEn  <- c(pBEn, pwrN)
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
             method="scABE", regulator=reg
             )
  class(ret) <- c("pwrA", class(ret))
  return(ret)

}
