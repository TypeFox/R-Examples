nptestind <- function(formula, data, subset,
                      rts = c("C", "NI", "V"), base = c("output", "input"),
                      reps = 999, alpha = 0.05,
                      print.level = 1, dots = TRUE){
 
 if (alpha < 0 | alpha > 99.99) {
  stop("'alpha' must be between 0 and 1 inclusive", call. = FALSE)
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
                       base = base, rts = rts, print.level = print.level,
                       type = "DF", winw = winw)
 
 # get the data in matrices
 
 Y  <- YX$y
 X  <- YX$x
 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 rt <- YX$myrts
	ba <- YX$mybase
	esample <- YX$esample

	# original Farrell measures

	te <- .teRad(t(Y),t(X),M,N,K,t(Y),t(X),K,rt,ba,1,print.level=0)
	
	# redefine if some Farrell measures are not computed
	
	te.good <- !is.na(te)
	K  <- sum(te.good)
	if(K == 0){
  stop("Could not compute measure of technical efficiency for a single data point")
	}
	te <- te[te.good]
	Y  <- Y[te.good, , drop = FALSE]
	X  <- X[te.good, , drop = FALSE]
	esample[!te.good] <- FALSE
	
	te <- ifelse(abs(te-1) < .Machine$double.eps, 1, te)
	
	# Begin Test
	
	# step 1
	
	terfl <- c(te, 2-te)
	mybw  <- bw.SJ(terfl, method = c("dpi"))
	scVarHom <- ( 1 + mybw^2 / var(terfl) )^(-1/2)
	
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
	
	# step 2
	t4n = .t4n(w = Z1, d = te, FALSE)
	# print(t4n)
	# print(Z1)
	
	# step 3
	# done
	
 # begin bootstrap
	
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

	if(print.level >= 1){
	 if ( ba == 1){
	  if(N == 1){
	   inps <- "input"
	  }
	  else {
	   inps <- "mix of inputs"
	  }
# 	  cat(" Test: Ho: input-based measure of technical efficiency and\n", sep = "")
# 	  cat("           ",inps," are independent \n\n", sep = "")
	  cat("", paste("",rep("_", (winw-10)/1),"", sep = ""), sep = "")
	  cat("\n          Test\n\n", sep = "")
	  mymesage <- paste("Ho: input-based measure of technical efficiency and ",inps," are independent\n", sep = "")
	 }
	 else {
	  if(M == 1){
	   outs <- "output"
	  }
	  else {
	   outs <- "mix of outputs"
	  }
# 	  cat(" Test: Ho: output-based measure of technical efficiency and\n", sep = "")
# 	  cat("           ",outs," are independent \n\n", sep = "")
	  mymesage <- paste("\nHo: input-based measure of technical efficiency and ",outs," are independent\n", sep = "")
	 }
	 cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
	}
	 
 
	t4nb <- numeric(reps)
	boot.type <- paste(" Bootstrapping test statistic T4n (",reps," replications)\n", sep = "")
	mychar <- "."
	
 for(b in seq_len(reps)){
  if (dots){
   if (winw != 0){
    if(b == 1) .dots(0, boot.type, width = winw)
   }
  }
  # step 4
  teb <- .smplHomTE(terfl,Kr=K,mybw,scVarHom,ba)
  # step 5
  sample_w <- floor(K * runif(K) + 1)
  wb <- Z1[sample_w,]
  # step 6
  t4nb[b] <- .t4n(w = wb, d = teb)
  # dots
  if (dots){
   if (winw != 0){
    .dots(b, width = winw, character = mychar)
   }
  }
 } # step 7
	# cat("\n")
	
	# step 8 

	pval <- mean(t4nb >= t4n, na.rm = TRUE)
	
	if(print.level >= 1){
	 if(ba ==1){
#    cat(" p-value of the Ho that input-based measure of technical efficiency and\n", sep = "")
#    cat("                        ",inps," are independent = ",formatC(pval, digits = 4, format = "f"),"\n", sep = "")
   mymesage <- paste("\np-value of the Ho that input-based measure of technical efficiency and ",inps," are independent = ",formatC(pval, digits = 4, format = "f"),":", sep = "")
	 }
	 else {
#    cat(" p-value of the Ho that output-based measure of technical efficiency and\n", sep = "")
#    cat("                        ",outs," are independent = ",formatC(pval, digits = 4, format = "f"),"\n", sep = "")
   mymesage <- paste("\np-value of the Ho that output-based measure of technical efficiency and ",outs," are independent = ",formatC(pval, digits = 4, format = "f"),":", sep = "")
	 }
	 cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
	
	 mymesage <- paste("\n",ifelse(pval <= alpha, "Heterogeneous", "Homogeneous")," bootstrap ",ifelse(pval <= alpha, "should", "can")," be used when performing ",YX$base.string,"-based technical efficiency measurement under assumption of ",YX$rts.string," technology\n", sep = "")
	 cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
# 	if(pval <= alpha){
# 	 mymesage <- paste("\nHeterogeneous bootstrap \n", sep = "")
# 	 cat("\n\n");
# 	 cat("   d and w are not independent ","\n");
# 	 cat("   go ahead with SW 2000 heterogeneous bootstrap ","\n\n\n");
# 	 decision = "d and w are not independent ==> SW 2000 heterogeneous bootstrap"
# 	} else {
# 	 cat("\n\n");
# 	 cat("   d and w are independent ","\n");
# 	 cat("   go ahead with SW 1998 homogeneous bootstrap ","\n\n\n");
# 	 decision = "d and w are independent ==> SW 1998 homogeneous bootstrap"
# 	}
	}
	tymch <- list(K = K, M = M, N = N, reps = reps, alpha = alpha,
	              t4n = t4n, pval = pval)
	class(tymch) <- "npsf"
	return(tymch)

}



#