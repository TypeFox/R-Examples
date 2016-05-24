teradialbc <- function(formula, data, subset,
                       ref = NULL, data.ref = NULL, subset.ref = NULL,
                       rts = c("C", "NI", "V"), base = c("output", "input"),
                       homogeneous = TRUE, smoothed = TRUE, kappa = NULL,
                       reps = 999, level = 95,
                       core.count = 1, cl.type = c("SOCK", "MPI"),
                       print.level = 1, dots = TRUE){
 if( !is.null(ref) & is.null(data.ref) ){
  warning("If you use variable names in 'ref', 'data.ref' is required", call. = FALSE)
 }

 if( is.null(ref) & !is.null(data.ref) ){
  stop("If you use 'data.ref', 'ref' is required", call. = FALSE)
 }
 
 if (level < 10 | level > 99.99) {
  stop("'level' must be between 10 and 99.99 inclusive")
 }
 
 if(!is.null(kappa)){
  if (kappa <= 0.5 | kappa >= 1) {
   stop("'kappa' must be between 0 and 1")
  }
 }
 else {
  if(!smoothed){
   stop("'kappa' must be provided for subsampling bootstrap")
  }
 }

 if (reps < 100) {
  stop("'reps' must be at least 100")
 }
 
 if (reps < 200) {
  warning(" Statistical inference may be unreliable \n          for small number of bootstrap replications", call. = FALSE, immediate. = TRUE)
  warning(" Statistical inference may be unreliable for small number of bootstrap replications\n", call. = FALSE, immediate. = FALSE)
 }
 
 if (reps > 2000) {
  warning(" Unnecessary too many bootstrap replications; \n          consider setting 'reps' smaller than 2000", call. = FALSE, immediate. = TRUE)
  warning(" Unnecessary too many bootstrap replications; consider setting 'reps' smaller than 2000\n", call. = FALSE, immediate. = FALSE)
 }
 
 # begin require for parallel computing
 if (!is.numeric(core.count) || floor(core.count) != core.count || 
     core.count < 1){
  stop("'core.count' must be a positive integer", call. = FALSE)
 }
 
 if (core.count > 1){
  if (!(is.character(cl.type)) || !(cl.type %in% c("SOCK", "MPI"))){
   stop("invalid cluster type in 'cl.type'; 'cl.type' must be \"SOCK\" or \"MPI\"", call. = FALSE)
  }
  if (!("snowFT" %in% rownames(installed.packages()))){
   # mymessage <- unlist(strsplit("Package 'snowFT' required for parallel computing is not installed; type -install.packages(''snowFT'')- to install it", split = " "))
   stop("Package 'snowFT' required for parallel computing is not installed; \nuse -install.packages(",paste(dQuote("snowFT")), ")- to install it\n")
  }else{
   # require("snowFT")
   requireNamespace("snowFT", quietly = TRUE)
  }
  if (cl.type == "MPI"){
   if (!("Rmpi" %in% rownames(installed.packages()))){
    stop("Package 'Rmpi' required for parallel computing with option 'MPI' is not installed; \nuse -install.packages(",paste(dQuote("Rmpi")), ")- to install it\n")
   }else{
    # require("Rmpi")
    requireNamespace("Rmpi", quietly = TRUE)
   }
  }
 }
 # end require for parallel computing

 winw <- getOption("width")
 if (winw > 100+5){
  winw <- 100
 }
 else if (winw > 90+5) {
  winw <- 90
 }
 else if (winw > 80+5) {
  winw <- 80
 }
 else if (winw > 70+5) {
  winw <- 70
 }
 else if (winw > 60+5) {
  winw <- 60
 }
 else if (winw > 50+5) {
  winw <- 50
 }
 else {
  winw <- 0
 }
 
 # get the data in matrices
  
 YX <- .prepareYX(formula = formula, data = data, subset = subset, rts = rts,
                  base = base, ref = ref,	data.ref = data.ref, subset.ref = subset.ref,
                  print.level = print.level, type = "DF", winw = winw)
 Y  <- YX$y
 X  <- YX$x
 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 Yr <- YX$y.ref
	Xr <- YX$x.ref
	Kr <- nrow(Yr)
 rt <- YX$myrts
	ba <- YX$mybase
	esample <- YX$esample
	
	# original Farrell measures
	
	te0 <- .teRad(t(Y),t(X),M,N,K,t(Yr),t(Xr),Kr,rt,ba,1,print.level=print.level)
	
	# redefine if some Farrell measures are not computed
	
	te.good <- !is.na(te0)
	K  <- sum(te.good)
	if(K == 0){
  stop("Could not compute measure of technical efficiency for a single data point")
	}
	te <- te0[te.good]
	Y  <- Y[te.good, , drop = FALSE]
	X  <- X[te.good, , drop = FALSE]
	esample[!te.good] <- FALSE
	
	
	# begin preparation for bootstrap
	
 if( smoothed ){
  
  # Common for homogeneous and heterogeneous
  
  # Calculate Farrell measures for reference data points: teRef
  
  if( is.null(ref) ){
   teRef <- te
  }
  else {
   teRef0 <- .teRad(t(Yr),t(Xr),M,N,Kr,t(Yr),t(Xr),Kr,rt,ba,
                    1,print.level=0)
   
   # redefine if some Farrell measures are not computed
   
   teRef.good <- !is.na(teRef0)
   teRef <- teRef0[teRef.good]
   Yr    <- Yr[teRef.good, , drop = FALSE]
   Xr    <- Xr[teRef.good, , drop = FALSE]
   Kr    <- sum(teRef.good)
  }
  
  terfl <- c(teRef, 2-teRef)
  mybw  <- bw.SJ(terfl, method = c("dpi"))
  
  msub <- NULL
	 
  if ( homogeneous ){
   # begin homogeneous
   # print(1)
   scVarHom <- ( 1 + mybw^2 / var(terfl) )^(-1/2)
   # end homogeneous
  }
  else {
   # begin heterogeneous
   # print(2)
   if( ba == 2 ) {
    # output
    if( M == 1 ) {
     Z1 <- cbind(Yr, Xr)
    }
    else {
     nu <- atan( Yr[ , -1, drop = FALSE] / Yr[ , 1] )
     pot.zeros <- Yr[ , 1] == 0
     if( any(pot.zeros) ){
      nu[pot.zeros, ] <- matrix(0, nrow = sum(pot.zeros), ncol = M-1)
     }
     Z1 <- cbind( nu, Xr )
    }
   }
   else {
    # input
    if( N == 1 ) {
     Z1 <- cbind(Yr, Xr)
    }
    else {
     nu <- atan( Xr[ , -1, drop = FALSE] / Xr[ , 1] )
     pot.zeros <- Xr[ , 1] == 0
     if( any(pot.zeros) ){
      nu[pot.zeros, ] <- matrix(0, nrow = sum(pot.zeros), ncol = N-1)
     }
     Z1 <- cbind( Yr, nu )
    }
   }
   
   bws <- apply(Z1, MARGIN = 2, bw.SJ, method = c("dpi"))
   bws <- matrix( c(bws, mybw), nrow = 1)
   scVarHet <- (1 + bws^2)^(-1/2)
   
   Z  <- cbind( Z1, teRef )
   Zr <- cbind( Z1, 2-teRef )
   L1 <- chol( cov(Z) )
			L2 <- chol( cov(Zr) )
   Zt <- rbind( Z,Zr )
   nZ <- ncol(Zt)
   
   onesN <- matrix(1, nrow = Kr, ncol = 1)
   M1    <- diag(Kr) - matrix(1/Kr, nrow = Kr, ncol = Kr)
   cons1 <- kronecker (onesN, scVarHet)
   cons2 <- kronecker (onesN, bws)
   # end heterogeneous
  }
 }
 else {
  # begin subsampling
  # print(3)
  msub  = round(Kr^kappa)
  # end subsampling
 }
	
 # end preparation for bootstrap	
	
 # begin bootstrap
	
	# printing
	if(print.level >= 1){
	 mymesage <- paste("\nBootstrapping reference set formed by ",Kr," ",ifelse(is.null(ref), "reference data", "data")," point(s) and computing radial (Debreu-Farrell) ",YX$base.string,"-based measures of technical efficiency under assumption of ",YX$rts.string," technology for each of ",K," data point(s) realtive to the bootstrapped reference set\n", sep = "")
	 cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
	 if( smoothed ){
   if ( homogeneous ) {
    # cat("\n Smoothed homogeneous bootstrap (",reps," replications)\n\n", sep = "")
    boot.type <- paste("Smoothed homogeneous bootstrap (",reps," replications)\n", sep = "")
   }
   else {
    # cat("\n Smoothed heterogeneous bootstrap (",reps," replications)\n\n", sep = "")
    boot.type <- paste("Smoothed heterogeneous bootstrap (",reps," replications)\n", sep = "")
   }
	 }
	 else {
	  # cat("\n Subsampling bootstrap (",reps," replications)\n\n", sep = "")
	  boot.type <- paste("Subsampling bootstrap (",reps," replications)\n", sep = "")
	 }
	}

	teboot <- matrix(nrow = reps, ncol = K)
	realnrep <- numeric(reps)
	
	# begin parallel computing
	if (core.count > 1){
	 if( smoothed ){
	  if ( homogeneous ) {
	   run.fun <- .run.boots.hom.bc
	   addInput <- list(Y = Y, X = X, Yr = Yr, Xr = Xr, rt = rt, 
	                    ba = ba, teRef = teRef, terfl = terfl, 
	                    mybw = mybw, scVarHom = scVarHom)
	  }
	  else{
	   run.fun <- .run.boots.het.bc
	   addInput <- list(Y = Y, X = X, Yr = Yr, Xr = Xr, rt = rt, 
	                    Zt = Zt, L1 = L1, L2 = L2, M1 = M1, 
	                    cons1 = cons1, cons2 = cons2, onesN = onesN, ba = ba)
	  }
	 }
	 else {
	  run.fun <- .run.boots.subs.bc
	  addInput <- list(Y = Y, X = X, Yr = Yr, Xr = Xr, rt = rt, 
	                   ba = ba, msub = msub)
	 }
	 inputs <- as.list(rep(0, reps))
	 if (dots){
	  if (winw != 0){
	   .dots(0, boot.type, width = winw)
	  }
	  outputs <- snowFT::performParallel(count = core.count, x = inputs, fun = run.fun, 
	                             cltype = cl.type, printfun = .run.dots, 
	                             printargs = list(width = winw), printrepl = 1, 
	                             ... = addInput)
	 }
	 else {
	  outputs <- snowFT::performParallel(count = core.count, x = inputs, fun = run.fun, 
	                             cltype = cl.type, 
	                             ... = addInput)
	 }
	 tymch3 <- matrix(unlist(outputs), nrow = reps, byrow = TRUE)
	 teboot <- tymch3[,1:K]
	 if ( !homogeneous ) {
	  realnrep <- tymch3[,K + 1]
	 }
	}
	# end parallel computing
	else {
	 for(b in seq_len(reps)){
	  if (print.level >= 1){
	   if (winw != 0){
	    if(b == 1) .dots(0, boot.type, width = winw)
	   }
	  }
	  if( smoothed ){
	   if ( homogeneous ) {
	    # begin homogeneous
	    # print(1)
	    # teB <- numeric(K)
	    if (ba == 2) {
	     # step 1: sampling
	     Yrb <- .smplHom(teRef,terfl,Kr,mybw,scVarHom,Yr,ba)
	     # step 2: applying DEA
	     teB <- .teRad(t(Y),t(X),M,N,K,t(Yrb),t(Xr),Kr,rt,ba,
	                   0,print.level=0)
	    }
	    else {
	     # step 1: sampling
	     Xrb <- .smplHom(teRef,terfl,Kr,mybw,scVarHom,Xr,ba)
	     # step 2: applying DEA
	     teB <- .teRad(t(Y),t(X),M,N,K,t(Yr),t(Xrb),Kr,rt,ba,
	                   0,print.level=0)
	    }
	    mychar <- "."
	    # end homogeneous
	   }
	   else {
	    # begin heterogeneous
	    # print(2)
	    # step 1: sampling
	    smplHet <- .smplHet(Zt,L1,L2,M1,cons1,cons2,onesN,Kr,nZ,ba,M,N)
	    # step 2: get efficiency under H0
	    teRefB  <- .teRad(t(smplHet$Yrb),t(smplHet$Xrb),M,N,smplHet$Krb,
	                      t(Yr),t(Xr),Kr,rt,ba,
	                      0,print.level=0)
	    # step 3: get efficient xRef or yRef
	    if (ba == 1) {
	     Yrb <- smplHet$Yrb
	     Xrb <- smplHet$Xrb * teRefB / smplHet$teb
	    }
	    else {
	     Yrb <- smplHet$Yrb * teRefB / smplHet$teb
	     Xrb <- smplHet$Xrb
	    }
	    # handle infeasible
	    teRefB.good <- !is.na(teRefB)
	    # modify the existing matrices by tossing the missing values
	    Yrb <- Yrb[teRefB.good, , drop = FALSE]
	    Xrb <- Xrb[teRefB.good, , drop = FALSE]
	    # new number of obs in reference
	    KrB <- sum(teRefB.good)
	    # step 4: get the bootstrapped TE
	    teB <- .teRad(t(Y),t(X),M,N,K,t(Yrb),t(Xrb),KrB,rt,ba,
	                  0,print.level=0)
	    # print(c(KrB,Kr,KrB/Kr))
	    if (KrB/Kr < 0.70){
	     mychar <- "?"
	    }
	    else if (KrB/Kr < 0.85){
	     mychar <- "x"
	    }
	    else {
	     mychar <- "."
	    }
	    realnrep[b] <- KrB
	    # end heterogeneous
	   }
	   # end smoothed
	  }
	  else {
	   # begin subsampling
	   # print(3)
	   # step 1: sub-sampling
	   newSample <- sample( seq_len(Kr), msub, replace = TRUE)
	   ystar <- Yr[newSample, , drop = FALSE]
	   xstar <- Xr[newSample, , drop = FALSE]
	   # step 2: applying DEA
	   teB <- .teRad(t(Y),t(X),M,N,K,t(ystar),t(xstar),msub,rt,ba,
	                 0,print.level=0)
	   # print(c(msub, Kr))
# 	   if (msub/Kr < 0.80){
# 	    mychar <- "?"
# 	   }
# 	   else if (msub/Kr < 0.95){
# 	    mychar <- "@"
# 	   }
# 	   else if (msub/Kr < 1){
# 	    mychar <- "x"
# 	   } 
# 	   else {
	    mychar <- "."
	   # }
	   # end subsampling
	  }
	  teboot[b,] <- teB
	  # print(teB[1:10])
	  
	  # dots
	  
	  if (dots){
	   if (winw != 0){
	    .dots(b, width = winw, character = mychar)
	   }
	  }
	 }
	}
	if(reps %% winw != 0) cat("\n")
	
	# Working with large BOOT matrix
	# Step 3: bias correction and CIs
	
	bci <- .biasAndCI(te, teboot, msub = msub, K, Kr, M, N,
	                 level, smoothed, forceLargerOne = FALSE)

	# Check if any bias corrected measure is negative;
 #	If yes, do inference in terms of Shephard
 #	distance functions (must be input oriented)
	if( min(bci$tebc, na.rm = TRUE) < 0 ){
  bci <- .biasAndCI(te, teboot, msub = msub, K, Kr, M, N,
                   level, smoothed, forceLargerOne = TRUE)
  te <- 1/te
  if(print.level >= 1){
   warning(" One or more bias-corrected Farrell ",YX$base.string,"-based measures of technical efficiency is negative.  Analysis will be done in terms of Shephard distance function, a reciprocal of the Farrell measure. \n", call. = FALSE, immediate. = FALSE)
  }
	}
	
	# Check if bias squared over var is larger for all
	if( min(bci$BovV, na.rm = TRUE) < 1){
  warning(" For one or more data points statistic [(1/3)*(bias)^2/var] is smaller than 1; bias-correction should not be used.\n", call. = FALSE, immediate. = FALSE)
	}
	
	if ( !homogeneous & smoothed ){
	 tymch <- list(K = K, M = M, N = N, reps = reps, level = level,
	               te = te, tebc = bci$tebc, biasboot = bci$bias, varboot = bci$vari,
	               biassqvar = bci$BovV, realreps = bci$reaB,
	               telow = bci$LL, teupp = bci$UU, teboot = teboot,
	               samples.b = realnrep, esample = esample)
	}
	else {
	 tymch <- list(K = K, M = M, N = N, reps = reps, level = level,
	               te = te, tebc = bci$tebc, biasboot = bci$bias, varboot = bci$vari,
	               biassqvar = bci$BovV, realreps = bci$reaB,
	               telow = bci$LL, teupp = bci$UU, teboot = teboot, esample = esample)
	}
	
	class(tymch) <- "npsf"
 return(tymch)
}



#