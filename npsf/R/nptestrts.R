nptestrts <- function(formula, data, subset,
                      base = c("output", "input"),
                      homogeneous = TRUE, test.two = TRUE,
                      reps = 999, alpha = 0.05,
                      core.count = 1, cl.type = c("SOCK", "MPI"),
                      print.level = 1, dots = TRUE){
 
 if (alpha < 0 | alpha > 99.99) {
  stop("'alpha' must be between 0 and 1 inclusive", call. = FALSE)
 }
 
 if (reps < 100) {
  stop("'reps' must be at least 100")
 }
 
 if (reps < 200) {
  warning(" Statistical inference may be unreliable \n          for small number of bootstrap replications", call. = FALSE, immediate. = TRUE)
  warning(" Statistical inference may be unreliable for small number of bootstrap replications", call. = FALSE, immediate. = FALSE)
 }
 
 if (reps > 2000) {
  warning(" Unnecessary too many bootstrap replications; \n          consider setting 'reps' smaller than 2000", call. = FALSE, immediate. = TRUE)
  warning(" Unnecessary too many bootstrap replications; consider setting 'reps' smaller than 2000", call. = FALSE, immediate. = FALSE)
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
 
 YX <- .prepareYXnoRef(formula = formula, data = data, subset = subset,
                       base = base, print.level = print.level,
                       type = "DF", winw = winw, rts.subst = "CRS, NiRS, and VRS")
 Y  <- YX$y
 X  <- YX$x
 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 rt <- YX$myrts
	ba <- YX$mybase
	esample <- YX$esample

	# original Farrell measures

	# CRS
	teCrs <- .teRad(t(Y),t(X),M,N,K,t(Y),t(X),K,3,ba,1,print.level=0)

	# NIRS
	teNrs <- .teRad(t(Y),t(X),M,N,K,t(Y),t(X),K,2,ba,1,print.level=0)

	# VRS
	teVrs <- .teRad(t(Y),t(X),M,N,K,t(Y),t(X),K,1,ba,1,print.level=0)
	
	# Check for Missing Values
	teSum = teCrs + teNrs + teVrs
	
	# redefine if some Farrell measures are not computed
	
	te.good <- !is.na(teSum)
	K  <- sum(te.good)
	if(K == 0){
  stop("Could not compute measure of technical efficiency for a single data point")
	}
	teCrs <- teCrs[te.good]
	teNrs <- teNrs[te.good]
	teVrs <- teVrs[te.good]
	Y <- Y[te.good, , drop = FALSE]
	X <- X[te.good, , drop = FALSE]
	esample[!te.good] <- FALSE
	
	seCrsMean <- mean(teCrs) / mean(teVrs)
	seCrs     <- teCrs /teVrs
	seNrsMean <- mean(teNrs) / mean(teVrs)
	seNrs     <- teNrs /teVrs
	
	# Begin Test #1: CRS vs VRS
	
	if(print.level >= 1){
	 cat("", paste("",rep("_", (winw-10)/1),"", sep = ""), sep = "")
	 cat("\n          Test #1\n\n", sep = "")
	 cat(" Ho: mean(F_i^CRS)/mean(F_i^VRS) = 1\n", sep = "")
	 cat("   and\n", sep = "")
	 cat(" Ho: F_i^CRS/F_i^VRS = 1 for each of ",K," data ", ngettext(K, "point", "point(s)"), "\n", sep = "")
	}
	
	# TE:	under Ho, RTS is CRS
	teRef <- teCrs
	
	# begin preparation for bootstrap
	
	terfl <- c(teRef, 2-teRef)
	mybw  <- bw.SJ(terfl, method = c("dpi"))
	
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
	   Z1 <- cbind(Y, X)
	  }
	  else {
	   nu <- atan( Y[ , -1, drop = FALSE] / Y[ , 1] )
	   pot.zeros <- Y[ , 1] == 0
	   if( any(pot.zeros) ){
	    nu[pot.zeros, ] <- matrix(0, nrow = sum(pot.zeros), ncol = M-1)
	   }
	   Z1 <- cbind( nu, X )
	  }
	 }
	 else {
	  # input
	  if( N == 1 ) {
	   Z1 <- cbind(Y, X)
	  }
	  else {
	   nu <- atan( X[ , -1, drop = FALSE] / X[ , 1] )
	   pot.zeros <- X[ , 1] == 0
	   if( any(pot.zeros) ){
	    nu[pot.zeros, ] <- matrix(0, nrow = sum(pot.zeros), ncol = N-1)
	   }
	   Z1 <- cbind( Y, nu )
	  }
	 }
	 
	 bwsZ <- apply(Z1, MARGIN = 2, bw.SJ, method = c("dpi"))
	 bws <- matrix( c(bwsZ, mybw), nrow = 1)
	 scVarHet <- (1 + bws^2)^(-1/2)
	 
	 Z  <- cbind( Z1, teRef )
	 Zr <- cbind( Z1, 2-teRef )
	 L1 <- chol( cov(Z) )
	 L2 <- chol( cov(Zr) )
	 Zt <- rbind( Z,Zr )
	 nZ <- ncol(Zt)
	 
	 onesN <- matrix(1, nrow = K, ncol = 1)
	 M1    <- diag(K) - matrix(1/K, nrow = K, ncol = K)
	 cons1 <- kronecker (onesN, scVarHet)
	 cons2 <- kronecker (onesN, bws)
	 # end heterogeneous
	}

	# end preparation for bootstrap	
	
 # begin bootstrap
	
	# printing
	if(print.level >= 1){
	 mymesage <- paste("\nBootstrapping reference set formed by ",K," data point(s) and computing radial (Debreu-Farrell) ",YX$base.string,"-based measures of technical efficiency under assumptions of CRS and VRS technology for each of ",K," data point(s) realtive to the bootstrapped reference set\n", sep = "")
	 cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
# 	 cat(" Bootstrapping reference set formed by ",K," data points and computing \n", sep = "")
# 	 cat(" radial (Debreu-Farrell) ",YX$base.string,"-based measures of technical efficiency  \n", sep = "")
#   cat(" under assumptions of CRS and VRS technology for each of ",K," data points \n", sep = "")
#   cat(" realtive to the bootstrapped reference set \n", sep = "")
  if ( homogeneous ) {
    # cat("\n Smoothed homogeneous bootstrap (",reps," replications)\n\n", sep = "")
    boot.type <- paste("Smoothed homogeneous bootstrap (",reps," replications)\n", sep = "")
   }
   else {
    # cat("\n Smoothed heterogeneous bootstrap (",reps," replications)\n\n", sep = "")
    boot.type <- paste("Smoothed heterogeneous bootstrap (",reps," replications)\n", sep = "")
   }
	}
	
	
	te1boot <- te2boot <- matrix(nrow = reps, ncol = K)
	# te1booti <- te2booti <- numeric(K)
	realnrep <- numeric(reps)
	
	# begin parallel computing
	if (core.count > 1){
	 if ( homogeneous ) {
	  run.fun <- .run.boots.hom.rts
   addInput <- list(Y = Y, X = X, rtsHo = 3, ba = ba, teRef = teRef, terfl = terfl,
                    mybw = mybw, scVarHom = scVarHom)
  }else{
   run.fun <- .run.boots.het.rts
   addInput <- list(Y = Y, X = X, Yr = Y, Xr = X, 
                    rtsHo = 3, Zt = Zt, L1 = L1, L2 = L2, M1 = M1,
                    cons1 = cons1, cons2 = cons2, onesN = onesN, ba = ba)
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
	 }else{
	  outputs <- snowFT::performParallel(count = core.count, x = inputs, fun = run.fun, 
	                             cltype = cl.type, 
	                             ... = addInput)
  }
	 tymch3 <- matrix(unlist(outputs), nrow = reps, byrow = TRUE)
	 te1boot <- tymch3[,1:K]
	 te2boot <- tymch3[,(K + 1):(2 * K)]
	 if ( !homogeneous ) {
	  realnrep <- tymch3[,2 * K + 1]
	 }
	}
	# end parallel computing
	else {
  for(b in seq_len(reps)){
   if (dots){
    if (winw != 0){
     if(b == 1) .dots(0, boot.type, width = winw)
    }
   }
   for(ii in seq_len(K)){
     if ( homogeneous ) {
      # begin homogeneous
      # print(1)
      if (ba == 2) {
       # step 1: sampling
       Yb <- .smplHom(teRef,terfl,Kr=K,mybw,scVarHom,Y,ba)
       # step 2: applying DEA
       # CRS
       te1boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yb),t(X),K,3,ba,
                               0,print.level=0)
       # VRS
       te2boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yb),t(X),K,1,ba,
                               0,print.level=0)
      }
      else {
       # step 1: sampling
       Xb <- .smplHom(teRef,terfl,Kr=K,mybw,scVarHom,X,ba)
       # step 2: applying DEA
       # CRS
       te1boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Y),t(Xb),K,3,ba,
                               0,print.level=0)
       # VRS
       te2boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Y),t(Xb),K,1,ba,
                               0,print.level=0)
      }
      mychar <- "."
      # end homogeneous
     }
     else {
      # begin heterogeneous
      # print(2)
      # step 1: sampling
      smplHet <- .smplHet(Zt,L1,L2,M1,cons1,cons2,onesN,K,nZ,ba,M,N)
      # step 2: get efficiency under H0
      teRefB  <- .teRad(t(smplHet$Yrb),t(smplHet$Xrb),M,N,smplHet$Krb,
                       t(Y),t(X),K,rts=3,ba,
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
      # CRS
      te1boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yrb),t(Xrb),KrB,3,ba,
                     0,print.level=0)
      # VRS
      te2boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yrb),t(Xrb),KrB,1,ba,
                              0,print.level=0)
      
      if (KrB/K < 0.80){
       mychar <- "?"
      }
      else if (KrB/K < 0.95){
       mychar <- "@"
      }
      else if (KrB/K < 1){
       mychar <- "x"
      }
      else {
       mychar <- "."
      }
      realnrep[b] <- KrB
      # end heterogeneous
    }
    # end for i
   }
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
	
	# p-value calculation for Test #1
	pvalsTestOne <- .pvalsTestOne (seCrs, seCrsMean, te1boot, te2boot, ba)
	alphaL       <- 1 - (1-alpha)^(1/K)
	s.efficient  <- pvalsTestOne$plCRS > alphaL
	# s.efficient <- c(FALSE, rep(TRUE, K-1))
	n.sefficient <- sum(s.efficient, na.rm = TRUE)
	locComputed <- TRUE
	if(n.sefficient == 0){
  locComputed <- FALSE
  n.sefficient <- K
	}
	# print(s.efficient)
 # print(pvalsTestOne)
 # print(n.sefficient)
	s.inefficient <- ifelse(is.na(s.efficient), FALSE, !s.efficient)
	# print(s.inefficient)
	
	if(print.level >= 1){
	 mymesage <- paste("\np-value of the Ho that mean(F_i^CRS)/mean(F_i^VRS) = 1 (Ho that the global technology is CRS) = ",formatC(pvalsTestOne$pgCRS, digits = 4, format = "f"),":\n\nmean(hat{F_i^CRS})/mean(hat{F_i^VRS}) = ",formatC(seCrsMean, digits = 4, format = "f")," ",ifelse(pvalsTestOne$pgCRS>alpha, "is not", "is")," statistically greater than 1 at the ",alpha*100,"% significance level", sep = "")
	 # print((strsplit(mymesage, " ")))
	 cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
	 if(n.sefficient == K & locComputed){
	  cat("\nAll data points are scale efficient\n", sep = "")
	 }
	 if(!locComputed){
	  cat("\n", sep = "")
	  warning(" Statistical inference on individual scale efficiency is not \n          performed likely due to few bootstrap replications \n          (for all data points)", call. = FALSE, immediate. = TRUE)
	  warning(" Statistical inference on individual scale efficiency is not \n performed likely due to few bootstrap replications (for all data points)", call. = FALSE, immediate. = FALSE)
	 }
	 if(locComputed & min(pvalsTestOne$nonNa, na.rm = TRUE) < 100){
	  warning(" Statistical inference on individual scale efficiency is not \n          performed likely due to few bootstrap replications \n          (for some data points)", call. = FALSE, immediate. = TRUE)
	  warning(" Statistical inference on individual scale efficiency is not \n performed likely due to few bootstrap replications (for some data points)", call. = FALSE, immediate. = FALSE)
	 }
	}
	
	# Test #2: NIRS vs VRS
	
	if ( test.two ){
	 mychar <- "."
	 
	 n.sinefficient <- sum(s.inefficient)

	 if (pvalsTestOne$pgCRS > alpha) {
	  performGlobal <- FALSE
	 }
	 else {
	  performGlobal <- TRUE
	 }
	 
	 # do it only if
	 if(performGlobal || n.sinefficient > 0){
	  # Change the matrices, since now technology is NIRS under Ho; under Ho, RTS is NiRS
   teRef <- teNrs
	
   # begin preparation for bootstrap
	
   terfl <- c(teRef, 2-teRef)
   mybw  <- bw.SJ(terfl, method = c("dpi"))
	
   if ( homogeneous ){
    # begin homogeneous
    scVarHom <- ( 1 + mybw^2 / var(terfl) )^(-1/2)
    # end homogeneous
   }
   else {
    # begin heterogeneous
    bws <- matrix( c(bwsZ, mybw), nrow = 1)
    scVarHet <- (1 + bws^2)^(-1/2)
    Z  <- cbind( Z1, teRef )
    Zr <- cbind( Z1, 2-teRef )
    L1 <- chol( cov(Z) )
    L2 <- chol( cov(Zr) )
    Zt <- rbind( Z,Zr )
    cons1 <- kronecker (onesN, scVarHet)
    cons2 <- kronecker (onesN, bws)
    # end heterogeneous
   }
   # end preparation for bootstrap
	 }
	 
  # Begin bootstrap to perform Test #2: NIRS vs VRS

  if (pvalsTestOne$pgCRS > alpha){
   # CRS is not rejected so perform only individual tests
   if ( n.sinefficient > 0 ){
    if(print.level >= 1){
     cat("", paste("",rep("_", (winw-10)/1),"", sep = ""), sep = "")
     cat("\n          Test #2\n\n", sep = "")
     cat(" Ho: F_i^CRS/F_i^VRS = 1 for each of ",n.sinefficient," data ", ngettext(K, "point", "point(s)"), "\n", sep = "")
     mymesage <- paste("\nBootstrapping reference set formed by ",K," data point(s) and computing radial (Debreu-Farrell) ",YX$base.string,"-based measures of technical efficiency under assumptions of CRS and VRS technology for each of ",n.sinefficient," data point(s) realtive to the bootstrapped reference set\n", sep = "")
     cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
     
#      cat(" Bootstrapping reference set formed by ",K," data points and computing \n", sep = "")
#      cat(" radial (Debreu-Farrell) ",YX$base.string,"-based measures of technical efficiency  \n", sep = "")
#      cat(" under assumptions of NiRS and VRS technology for each of ",n.sinefficient," data points \n", sep = "")
#      cat(" realtive to the bootstrapped reference set \n", sep = "")
    }
    # begin case 1b
    Ysineff <- Y[s.inefficient, ,drop = FALSE]
    Xsineff <- X[s.inefficient, ,drop = FALSE]
    te1boot <- te2boot <- matrix(nrow = reps, ncol = n.sinefficient)

    # te1booti <- te2booti <- numeric(K)
    realnrep <- numeric(reps)
    
    # begin parallel computing
    if (core.count > 1){
     if ( homogeneous ) {
      run.fun <- .run.boots.hom.rts
      addInput <- list(Y = Ysineff, X = Xsineff, rtsHo = 2, ba = ba, 
                       teRef = teRef, terfl = terfl, 
                       mybw = mybw, scVarHom = scVarHom)
     }else{
      run.fun <- .run.boots.het.rts
      addInput <- list(Y = Ysineff, X = Xsineff, Yr = Y, Xr = X, 
                       rtsHo = 2, Zt = Zt, L1 = L1, L2 = L2, M1 = M1, 
                       cons1 = cons1, cons2 = cons2, onesN = onesN, ba = ba)
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
     }else{
      outputs <- snowFT::performParallel(count = core.count, x = inputs, fun = run.fun, 
                                 cltype = cl.type, 
                                 ... = addInput)
     }
     tymch3 <- matrix(unlist(outputs), nrow = reps, byrow = TRUE)
     te1boot <- tymch3[,1:n.sinefficient]
     te2boot <- tymch3[,(n.sinefficient + 1):(2 * n.sinefficient)]
     if ( !homogeneous ) {
      realnrep <- tymch3[,2 * n.sinefficient + 1]
     }
    }
    # end parallel computing
    else {
     for(b in seq_len(reps)){
      for(ii in seq_len(n.sinefficient)){
       if ( homogeneous ) {
        # begin homogeneous
        # print(1)
        if (ba == 2) {
         # step 1: sampling
         Yb <- .smplHom(teRef,terfl,Kr=K,mybw,scVarHom,Y,ba)
         # step 2: applying DEA
         # NiRS
         te1boot[b,ii] <- .teRad(Ysineff[ii,],Xsineff[ii,],M,N,1,t(Yb),t(X),K,2,ba,
                                 0,print.level=0)
         # VRS
         te2boot[b,ii] <- .teRad(Ysineff[ii,],Xsineff[ii,],M,N,1,t(Yb),t(X),K,1,ba,
                                 0,print.level=0)
        }
        else {
         # step 1: sampling
         Xb <- .smplHom(teRef,terfl,Kr=K,mybw,scVarHom,X,ba)
         # step 2: applying DEA
         # NiRS
         te1boot[b,ii] <- .teRad(Ysineff[ii,],Xsineff[ii,],M,N,1,t(Y),t(Xb),K,2,ba,
                                 0,print.level=0)
         # VRS
         te2boot[b,ii] <- .teRad(Ysineff[ii,],Xsineff[ii,],M,N,1,t(Y),t(Xb),K,1,ba,
                                 0,print.level=0)
        }
        # mychar <- "."
        # end homogeneous
       }
       else {
        # begin heterogeneous
        # print(2)
        # step 1: sampling
        smplHet <- .smplHet(Zt,L1,L2,M1,cons1,cons2,onesN,K,nZ,ba,M,N)
        # step 2: get efficiency under H0
        teRefB  <- .teRad(t(smplHet$Yrb),t(smplHet$Xrb),M,N,smplHet$Krb,
                          t(Y),t(X),K,rts=2,ba,
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
        # NiRS
        te1boot[b,ii] <- .teRad(Ysineff[ii,],Xsineff[ii,],M,N,1,t(Yrb),t(Xrb),KrB,2,ba,
                                0,print.level=0)
        # VRS
        te2boot[b,ii] <- .teRad(Ysineff[ii,],Xsineff[ii,],M,N,1,t(Yrb),t(Xrb),KrB,1,ba,
                                0,print.level=0)
        # end heterogeneous
       }
       # end for ii
      }
      # dots
      if (dots){
       if (winw != 0){
        if(b == 1) .dots(0, boot.type, width = winw)
        .dots(b, width = winw, character = mychar)
       }
      }
      # end for b
     }
    } # end with one core
    cat("\n")
   } #end of if ( n.sinefficient > 0 ){
  }
  else {
   # CRS is rejected, perform NiRS global and individual tests
   if(print.level >= 1){
    cat("", paste("",rep("_", (winw-10)/1),"", sep = ""), sep = "")
    cat("\n          Test #2\n\n", sep = "")
    cat(" Ho: mean(F_i^NiRS)/mean(F_i^VRS) = 1\n", sep = "")
    if (n.sinefficient > 0) {
     cat("   and\n", sep = "")
     cat(" Ho: F_i^NiRS/F_i^VRS = 1 for each of ",n.sinefficient," data ", ngettext(K, "point", "point(s)"), "\n", sep = "")
    }
    mymesage <- paste("\nBootstrapping reference set formed by ",K," data point(s) and computing radial (Debreu-Farrell) ",YX$base.string,"-based measures of technical efficiency under assumptions of CRS and VRS technology for each of ",K," data point(s) realtive to the bootstrapped reference set\n", sep = "")
    cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   }
   # begin case 2a
   te1boot <- te2boot <- matrix(nrow = reps, ncol = K)
   # te1booti <- te2booti <- numeric(K)
   realnrep <- numeric(reps)
   
   # begin parallel computing
   if (core.count > 1){
    if ( homogeneous ) {
     run.fun <- .run.boots.hom.rts
     addInput <- list(Y = Y, X = X, rtsHo = 2, ba = ba, teRef = teRef, terfl = terfl, 
                      mybw = mybw, scVarHom = scVarHom)
    }else{
     run.fun <- .run.boots.het.rts
     addInput <- list(Y = Y, X = X, Yr = Y, Xr = X,
                      rtsHo = 2, Zt = Zt, L1 = L1, L2 = L2, M1 = M1, 
                      cons1 = cons1, cons2 = cons2, onesN = onesN, ba = ba)
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
    }else{
     outputs <- snowFT::performParallel(count = core.count, x = inputs, fun = run.fun, 
                                cltype = cl.type, 
                                ... = addInput)
    }
    tymch3 <- matrix(unlist(outputs), nrow = reps, byrow = TRUE)
    te1boot <- tymch3[,1:K]
    te2boot <- tymch3[,(K + 1):(2 * K)]
    if ( !homogeneous ) {
     realnrep <- tymch3[,2 * K + 1]
    }
   }
   # end parallel computing
   else {
    for(b in seq_len(reps)){
     for(ii in seq_len(K)){
      if ( homogeneous ) {
       # begin homogeneous
       # print(1)
       if (ba == 2) {
        # step 1: sampling
        Yb <- .smplHom(teRef,terfl,Kr=K,mybw,scVarHom,Y,ba)
        # step 2: applying DEA
        # NiRS
        te1boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yb),t(X),K,2,ba,
                                0,print.level=0)
        # VRS
        te2boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yb),t(X),K,1,ba,
                                0,print.level=0)
       }
       else {
        # step 1: sampling
        Xb <- .smplHom(teRef,terfl,Kr=K,mybw,scVarHom,X,ba)
        # step 2: applying DEA
        # NiRS
        te1boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Y),t(Xb),K,2,ba,
                                0,print.level=0)
        # VRS
        te2boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Y),t(Xb),K,1,ba,
                                0,print.level=0)
       }
       # mychar <- "."
       # end homogeneous
      }
      else {
       # begin heterogeneous
       # print(2)
       # step 1: sampling
       smplHet <- .smplHet(Zt,L1,L2,M1,cons1,cons2,onesN,K,nZ,ba,M,N)
       # step 2: get efficiency under H0
       teRefB  <- .teRad(t(smplHet$Yrb),t(smplHet$Xrb),M,N,smplHet$Krb,
                         t(Y),t(X),K,rts=2,ba,
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
       # NiRS
       te1boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yrb),t(Xrb),KrB,2,ba,
                               0,print.level=0)
       # VRS
       te2boot[b,ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yrb),t(Xrb),KrB,1,ba,
                               0,print.level=0)
       # end heterogeneous
      }
      # end for ii
     }
     # dots
     if (dots){
      if (winw != 0){
       if(b == 1) .dots(0, boot.type, width = winw)
       .dots(b, width = winw, character = mychar)
      }
     }
    } # end for b
   } # end with one core
  } # end of CRS is rejected
	 if(reps %% winw != 0) cat("\n")
	 if ( n.sinefficient == 1){
	  te1boot <- matrix(te1boot, ncol = 1)
	  te2boot <- matrix(te2boot, ncol = 1)
	 }
# 	 print(te1boot)
# 	 print(te2boot)
# 	 print(class(te1boot))
# 	 print(class(te2boot))
	 
	 # p-value calculation for Test #2
	 if(performGlobal || n.sinefficient > 0){
	  pvalsTestTwo <- .pvalsTestTwo (seNrs, seNrsMean, te1boot, te2boot, ba, performGlobal, s.inefficient)
	  sineffdrs <-  pvalsTestTwo$pineffdrs > alphaL
	  nsineffdrs <- sum(sineffdrs)
	  if(performGlobal){
	   if (print.level >= 1){
	    mymesage <- paste("\np-value of the Ho that mean(F_i^NiRS)/mean(F_i^VRS) = 1 (Ho that the global technology is NiRS) = ",formatC(pvalsTestTwo$pgNRS, digits = 4, format = "f"),":\n\nmean(hat{F_i^NiRS})/mean(hat{F_i^VRS}) = ",formatC(seNrsMean, digits = 4, format = "f")," ",ifelse(pvalsTestTwo$pgNRS>alpha, "is not", "is")," statistically greater than 1 at the ",alpha*100,"% significance level", sep = "")
	    # print((strsplit(mymesage, " ")))
	    cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
	   }
	  }
	 }
	} # end of if ( test.two )
	
	tymch <- list(K = K, M = M, N = N, reps = reps, alpha = alpha,
	              teCrs = teCrs, teNrs = teNrs, teVrs = teVrs, 
	              sefficiency = seCrs, sefficiencyMean = seCrsMean, 
	              pGlobalCRS = pvalsTestOne$pgCRS )
	if (locComputed) {
  tymch$psefficient <- pvalsTestOne$plCRS
  tymch$sefficient <- s.efficient
  tymch$nsefficient <- n.sefficient
	}
 
	if(test.two){
	 if (performGlobal){
	  tymch$nrsOVERvrsMean <- seNrsMean
	  tymch$pGlobalNRS <- pvalsTestTwo$pgNRS
	 }
	 if (n.sinefficient > 0){
	  tymch$sineffdrs <- sineffdrs
		 tymch$pineffdrs <- pvalsTestTwo$pineffdrs
		 tymch$nrsOVERvrs <- seNrs
	 }
	}
	tymch$esample < esample
	class(tymch) <- "npsf"
 return(tymch)
}



#