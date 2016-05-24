##########################################################################
#
#        
#    Outlier Detection for Multiplicative High-throughput Data
#
#             by
#
#      Soo-Heang Eo, PhD and HyungJun Cho, PhD
#      Deparment of Statistics 
#      Korea University
#
#      Last updated: December 2014
#
##########################################################################

###################################################################################
#
# Main function  for OutlierDM 
#
###################################################################################

setClass("OutlierDM", representation(call = "language",
									 raw.data = "data.frame",
									 res = "data.frame",
									 x.pair = "list",
									 k = "numeric",
									 outlier = "matrix",
									 n.outliers = "integer",
									 quantreg = "character",
									 method = "character",
									 contrl.para = "list"
									 )
	)

odm <- function(x, k = 3, quantreg = c("linear", "nonlin", "constant", "nonpar"), 
	method = c("proj", "diff", "pair", "grubbs", "dixon", "iqr",  "siqr", "Zscore"), ...) {
## Input parameters
# x : Data vectors or matrix
# k : tuning parameter 
# quantreg : a type of quantile regression models used for the modeling fitting
# method : a type of outlier detection methods
# ... : minor control parameters including pair.cre, dist.mthd, Lower, Upper, trans, and centering.
# pair.cre : creterion for pairwise approach
# dist.mthd : median or mean
# Lower : Lower values
# Upper : Upper values
# trans : a parameter for logarithm transformation.
# centering : Logical parameter for centering. If TRUE, data are centered by column means.
# projection.type : projection type, naive, PCA, LPC.
## output 

  ##########
  #Preparation
	call <- match.call()
	quantreg <- match.arg(quantreg)
	method <- match.arg(method)

	contrl.para = odm.control(...)
	pair.cre <- contrl.para[[1]]
	dist.mthd <- contrl.para[[2]]
	Lower <- contrl.para[[3]]
	Upper <- contrl.para[[4]]
	trans <- contrl.para[[5]]
	centering <- contrl.para[[6]]
	projection.type <- contrl.para[[7]]
	lbda <- contrl.para[[8]]
	nonlin.method <- contrl.para[[9]]
	nonlin.SS <- contrl.para[[10]]
	nonlin.Frank <- contrl.para[[11]]	
	ncl <- contrl.para[[12]]
	alpha <- contrl.para[[13]]

	### FIXME: what type of parallel computing is optimal for our function? ##	
	# Set parallel computing in order to incease computing power

	# Start constructing data matrix
	if(!is.data.frame(x)) x <- as.data.frame(x) 

  n.obs <- ncol(x) 
  n.para <- nrow(x) 
  rownames(x) <- 1:nrow(x) 
	#colnames(x) <- paste("N", 1:n.obs, sep = "")
	raw.data <- x
 
	if(n.para < 30) stop('The number of parameters (peptides) is too few.')

	# transformation
	x <- eval(call(trans,x))

	# End constructing data
	cat("Please wait... \n")

	if(method == "Zscore"){
	
		# Calssicial standard deviation criteria
		# Step 1. Compute the standard deviation s_j for each peptide j, and then z-score z_ij=(y_ij-y_(j))/s_j, where y_(j) and s_j are the sample mean and standard deviation.
		# Step 2. For each peptide j, sample i is flagged as an outlier if z_ij < -k  or z_ij > k, where k = 2 or 3. 
		#colnames(x) <- paste("Rep", 1:n.obs, sep = "")

		SD = apply(x, 1, sd)
		LB = -1 * k
		UB = k

		zij = t(apply(x, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
		colnames(zij) <- paste("Z_", colnames(zij), sep="")

		# Outputs
		x <- cbind(x, zij, SD, LB, UB)
		outlier = abs(zij) > k

		i = which( apply(outlier, 1, function(x) sum(x) != 0 ))
		n.outliers = length(i)

		Outlier = rep(FALSE, n.para)
		if(length(i) > 0) Outlier[i] <- TRUE
		x <- cbind(Outlier, x)	

		new("OutlierDM", 
			call = call, 
			raw.data = raw.data, 
			res = x, 
			outlier = outlier,
			n.outliers=n.outliers, 
			method = method, 
			k = k, 
			contrl.para = contrl.para)

	} else if(method =="dixon") {
  	  # Dixon's range test for multiplicative experiments
	   # The algorithm is based on the R package, outliers.
		  fit <- apply( x, 1, function(x) outliers::dixon.test(x)$p.value)
		  x <- cbind(x, pvalue = fit)
		  
		  i <-  which(x$pvalue <= 0.05)
		  n.outliers <- length(i)
		  Outlier <- rep(FALSE, n.para)
		  if(length(i) > 0) Outlier[i] <- TRUE
		  x <- cbind(Outlier, x)
		  
		  new("OutlierDM", 
		  	call = call, 
		  	raw.data = raw.data, 
		  	res = x, 
		  	n.outliers = n.outliers, 
		  	method = method, 
		  	contrl.para = contrl.para)

	} else if(method == "grubbs"){
  	# Grubbs test for multiplicative experiments
		Gij = apply(x, 1, function(x, alpha) rgrubbs.test(x), alpha = alpha)
		rownames(Gij) <- paste("G", rep(1:n.obs, each = 2), sep = "")
		
		g.sel = as.logical(1:nrow(Gij) %% 2)
    outlier = t(Gij[!g.sel,])
    #outlier <- ifelse(outlier < 0.05, TRUE, FALSE)
    outlier <- ifelse(is.na(outlier), FALSE, TRUE)

		#Outputs
		pval = apply( x, 1, function(x) outliers::grubbs.test(x)$p.value)
  		
		x <- cbind(x, t(Gij[g.sel,]))
      
   	i = which( pval <= alpha ) 
		n.outliers = length(i)

		Outlier <- rep(FALSE, n.para)
		if(length(i) > 0) Outlier[i] <- TRUE
		x <- cbind(Outlier, x)

		new("OutlierDM", 
			call 		= call, 
			raw.data = raw.data, 
			res 		= x, 
			outlier = outlier,
			n.outliers = n.outliers, 
			method 		= method, 
			contrl.para = contrl.para)

	} else if(method %in% c("iqr", "siqr")){
	  
    # compute the first and third quatiles for each peptide j and then its interquantile range IQRj
	  q1.j = apply(x, 1, function(x) quantile(x, prob = 0.25))
	  q2.j = apply(x, 1, function(x) quantile(x, prob = 0.50))
	  q3.j = apply(x, 1, function(x) quantile(x, prob = 0.75))
	  
    # interquantile range criteria and semi-interquantile range criteria
    if(method == "iqr"){
      iqr.j = q3.j - q1.j
      # For each j, sample i is flagged as an outlier if yij  Q1j - k * IQRj  
      LB = q1.j - k * iqr.j
      UB = q3.j + k * iqr.j
    } else{
      # for SIQR criteria
      siqrl.j = q2.j - q1.j 
      siqru.j = q3.j - q2.j

      # For ach j, sample i is flagged as an outlier if yij < Q1j - 2k * SIQRj or yij > Q3j + 2k * SIQRj
      LB = q1.j - 2*k*siqrl.j
      UB = q3.j + 2*k*siqru.j
    }
    
		# outputs
		outlier = x < LB | x > UB
		x <- cbind(x, Q1 = q1.j, Q2 = q2.j, Q3 = q3.j, LB = LB, UB = UB)
		
    i = which(apply(outlier, 1, function(x) sum(x) != 0))
    n.outliers = length(i)

    Outlier = rep(FALSE, n.para)
    if(length(i) > 0 ) Outlier[i] <- TRUE
    x <- cbind(Outlier, x)
        
    new("OutlierDM", 
    	call = call, 
    	raw.data = raw.data, 
    	res = x, 
    	outlier = outlier,
    	n.outliers = n.outliers, 
    	method = method, 
    	k = k, 
    	contrl.para = contrl.para)

	} else if(method == "pair"){
		# Pairwise outlier detection for multiplicative experiments
		# The algorithm is the same as that of OutlierD package when the number of replicaates are 2.
		# Step 1. Apply OutlierD to all possible pairs of experiments.
		# Step 2. Declare peptide j as an outlier if it is declared as an outlier for at least one pair.
		
  	nonlin.SS = "Asym" # fixed for nonlinear modelling

    fit <- switch(quantreg,
			constant = quant.const(x, Lower, Upper, method),
			linear = quant.linear(x, Lower, Upper, method),
			nonlin = quant.nonlin(x, Lower, Upper, method, nonlin.method, nonlin.SS, nonlin.Frank),
			nonpar = quant.nonpar(x, Lower, Upper, method, lbda)
		)
  
		#Outputs
		x <- cbind(x, 
			A = fit$A, 
			M = fit$M, 
			Q1 = fit$Q1, 
			Q3 = fit$Q3, 
			LB = fit$Q1-k*(fit$Q3-fit$Q1), 
			UB = fit$Q3+k*(fit$Q3-fit$Q1))

    j <-  which((x[,"M"] < x[,"LB"])|(x[,"M"] > x[,"UB"]))	
    n.outliers <- length(j)
    Outlier <- rep(FALSE, n.para)
    if(length(j) > 0) Outlier[j] <- TRUE
    x <- cbind(Outlier, x)
    
    #cat(" Done. \n")

		new("OutlierDM", 
			call = call, 
			raw.data = raw.data, 
			res 	= x, 
			k = k, 
			outlier = as.matrix(Outlier),
			n.outliers = n.outliers, 
			quantreg = quantreg, 
			method = method, 
			contrl.para = contrl.para)

	} else if(method == "diff"){ 
		# Difference approach for multiplicative highthroughput data 
		# Step 1. Compute the difference M = y - median and the average A .
		# Step 2. Obtain the first and third quantile values, Q1(A) and Q3(A), on a MA plot using quantile regression approach.
		# Step 3. Calculate IQR(A) = Q3(A) - Q1(A)
		# Step 4. Construct the lower and upper fences, LB(A) = Q1(A) - kIQR(A) and UB(A) = Q3(A) + kIQR(A), where k is a tuning parameter.
		# Step 5. Declare the i-th replicate of peptide j as an outlier if it locates over the upper fence or under the lower fence.
	  fit <- switch(quantreg,
	         constant = quant.const.D(x, dist.mthd, Lower, Upper, n.obs, n.para),
	         linear      = quant.linear.D(x, dist.mthd, Lower, Upper, n.obs, n.para),
	         nonlin     = quant.nonlin.D(x, dist.mthd, Lower, Upper, n.obs, n.para),
	         nonpar   = quant.nonpar.D(x, dist.mthd, Lower, Upper, n.obs, n.para, lbda)
	         )
	  #Outputs
	  new.x <- cbind(x, fit$M, fit$A, fit$Q1, fit$Q3, 
	          fit$Q1-k*(fit$Q3-fit$Q1), fit$Q3+k*(fit$Q3-fit$Q1))
	  colnames(new.x) <- c(paste("Rep",1:n.obs, sep = ""),paste("M",1:n.obs, sep = ""),"A","Q1","Q3","LB","UB")
	  i <- (new.x[,(n.obs+1):(2*n.obs)] < new.x$LB)|(new.x[,(n.obs+1):(2*n.obs)] > new.x$UB)
	  i <- apply(i,1, any)
	  n.outliers <- sum(i,na.rm = TRUE)
	  Outlier <- rep(FALSE, n.para)
	  if(sum(i,na.rm=TRUE) > 0) Outlier[i] <- TRUE
	  new.x <- cbind(Outlier, new.x)
	  cat("Done. \n")
		new("OutlierDM", 
			call = call, 
			raw.data = raw.data, 
			res = new.x, 
			k = k, 
			n.outliers = n.outliers, 
			quantreg = quantreg, 
			method = method, 
			contrl.para = contrl.para)
  } else if( method == "proj" ){ 
		# Outlier Detection using projections for multiplicative experiments
		# Step 1. Shift the sample means to the origin.(in preparation step)
		# Step 2. Find the first PC vector v using PCA on the space of y.(projection.type option)
		# Step 3. Obtain the projection of a vector of each peptide j on v.
		# Step 4. Compute the length of the projection A and the length of the difference between a vector of peptide j and the projection M.
		# Step 5. Obtain the third quantile value Q3(A), on a MA plot using a quantile regression approach, and assume Q1(A) = -Q3(A).
		# Step 6. Calculate IQR.
		# Step 7. Construct the lower and upper fences.
		# Step 8. Declare peptide j as an outlier if it locates over the upper fence or under the lower fence.

		#Centering
		if(centering){
			xMeans <- colMeans(x)
			x <- x - rep(xMeans, each = n.para)
		}
		
		if(projection.type == "naive"){
			pt.naive <- 1000 * rep(1,n.obs)
			pt.tmp <- t(apply( x, 1, dist.two, type = "Pred", q1 = rep(0,n.obs), q2 = pt.naive))
		}	else if(projection.type %in% c("rPCA","PCA")){
			if(projection.type == "rPCA"){ 
				# Project-Pursuit PCA
				pca1 <- PCAproj(x, method = "sd", CalcMethod = "lincomb")$loadings
				wt.line <- abs(pca1[,1])
			}	else if(projection.type == "PCA"){ # Classical PCA
				pca1 <- prcomp(x, scale=TRUE)
				wt.line <- abs(pca1$rotation[,1])
			}
			pt.pca <- 1000 * wt.line
			pt.tmp <- t(apply( x, 1,dist.two, type = "Pred", q1 = rep(0, n.obs), q2 = pt.pca))
		}	else{
			stop("TYPE input does not match in a function.")
		} 

		new.x <- cbind(M = pt.tmp[,1], A = pt.tmp[,2])

		fit <- switch( quantreg,
			   constant	= quant.const(new.x, Lower, Upper, method),
			   linear		= quant.linear(new.x,Lower, Upper, method),
			   nonlin 	= quant.nonlin(new.x,Lower, Upper, method, nonlin.method, nonlin.SS, nonlin.Frank),
			   nonpar 	= quant.nonpar(new.x,Lower, Upper, method, lbda)
			   )

		#Outputs
		x <- cbind(x, A = fit$A, M = fit$M, Q1 = fit$Q1, Q3 = fit$Q3, 
		    LB = fit$Q1-k*(fit$Q3-fit$Q1), UB = fit$Q3+k*(fit$Q3-fit$Q1))

		outlier = x$M >= x$UB

		i =  which(outlier)
		n.outliers = length(i)
		Outlier = rep(FALSE, n.para)
		if(length(i) >0) Outlier[i] <- TRUE

		x <- cbind(Outlier, x)
		#cat(" Done. \n")

		new("OutlierDM", 
			call = call, 
			raw.data = raw.data, 
			res 	=	x, 
			k 		=	k, 
			outlier = as.matrix(outlier),
			n.outliers = n.outliers, 
			quantreg = quantreg, 
			method = method, 
			contrl.para = contrl.para)
  }
}
#END################################################################s
