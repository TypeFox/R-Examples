
scad_L2.svc <- function(lambda1 = 0.01, lambda2=0.01, x, y, a = 3.7, tol = 10^(-4), class.weights= NULL, maxIter=1000, verbose=TRUE){
# SCAD SVM classification
#
# Input:
#   xtrain : n-by-d data matrix to train (n chips/patients, d clones/genes, d>>n )
#   ytrain : column vector of target {-1, 1}'s (for n chips/patiens )
#   a : tuning parameter in scad function (default: 3.7 or whatever the paper
#uses)
#   lambda1 : tuning parameter in scad function (default : 0.01 or whatever the
#paper uses)
# 	 lambda2 : tuning parameter in L2 function (default : 0.01) 
#   tol: the cut-off value to be taken as 0
#
# By Natalia Becker (08.08.2008) :-) or better: (20.08.2008) 
# 

# SVM mit variable selectoion (clone selection) using scad penalty.
# NB ! xtrain = t(EP-matrix) !!!! 
#res<-scadsvc (x=xtrain, y=ytrain, a = 3.7, tol = 10^(-3))
# print(res)

	# checks
	if(nrow(x) != length(y)) stop("Wrong dimensions: nrow(x) should be equal to length(y) !")
	if (nlevels(as.factor(y)) !=2 ) stop(paste("We need 2 classes, currently have:",  paste(levels(as.factor(y)), collapse=", ")) ) 
	
	xtrain <- x
	ytrain <- as.numeric(as.character(y))
	
	# start with linear svm:
	require(e1071, quietly=TRUE)
	require(corpcor) # for pseudoinverse
	require(statmod) # for matrixmultplications 
	
	# change class labels to -1 and 1
	ytrain <- 2*as.numeric(as.factor(y))-3 

	linsvc <- svm(xtrain, factor(ytrain), kernel="linear", class.weights=class.weights, fitted =FALSE)
			
	# index of the support vectors
	index <- linsvc$index
	# type of svm
	type <- linsvc$type
	# w: coefficients times the training labels * SV
	w = apply(as.vector(linsvc$coefs) * linsvc$SV, 2, sum)  
	
	# rho = Bias: A scalar value representing the bias or threshold of the SVM
	#classifier, 
	# which is the negative intercept.
	b = linsvc$rho
	diff = 1
	ntrain = nrow(xtrain)
	d = ncol(xtrain)   
	xind = 1:d
	i<-1
	if (verbose) print("start iterations:")
	while (diff > tol ) {
		
		
		# should write i in the same position
		#if (i %% 1000 !=0 ){
			#cat("\b\b\b\b\b",i)
			#flush.console()
		#}
		
		
		if (!is.null(maxIter)){ 
			if  (i > maxIter) {
				if (verbose) print(paste("max. iteration number of", maxIter  ,"has been reached. Stop iterations "))
				nw <- w
			  nb <- b
			
				break
			}	
		}
			
		x1 = cbind(rep(1,nrow(xtrain)), xtrain)
		y1 = ytrain
		sgnres = as.vector(y1 - x1 %*% c(b, w))
		# important!!!! : sometimes a point is lying exactly on the hyperline --> sgnres = 0
		# produce errors --> move this randomly at the one or the other size.
		sgnres[sgnres== 0] <- sample(c(1,-1),1) *  10^(-100)  
		res = abs(sgnres)
		y0 = y1 / res
		
		# ###
		#D = 1/(2*ntrain) * diag(1/res)
		# ###
		# save as a vector D_vec
		D_vec = 1/(2*nrow(xtrain)) * (1/res)
		aw = abs(w)
		dp = lambda1*(aw<=lambda1)+(a*lambda1-aw)/(a-1) * (lambda1<aw &  aw <=a*lambda1)
		Q1_vec<-c(0, dp/aw)
		
		#Q1 = diag( c(0, dp/aw))
		#Q = t(x1) %*% D %*% x1 + Q1
		# ps_Q<-pseudoinverse(Q)
		
		P = (1/(2* ntrain)) * t(y1 + y0) %*% x1 
		# by Scad + L2 use P.ext instead of P; P.ext= P - Q2
		
		#	Q3<-matrix(c(0, 2*lambda2), nrow=1 );
		#	colnames(Q3)<- c("", colnames(xtrain)) 
		
		# NEW	 Q ext= Q + Q3 Q3 = diag( c(0, 2*lambda2))
		#inv_Q<-.find.inverse(U=t(x1),D_vec=D_vec,A_vec=c(0, dp/aw + 2*lambda2))
		# NEW
		#nwb = inv_Q %*% t(P)
		nwb = .calc.mult.inv_Q_mat2(U=t(x1),D_vec=D_vec,A_vec=c(0, dp/aw + 2*lambda2), mat2=t(P), n.thr=500)	
				
		
		# for one case is not able to calculate nb, nw: result: Inf, -Inf , a lot of NaN
		# then abort the loop, save the coff from previous iteration 
		if ( any(!is.finite(as.numeric(names(table(nwb))))) ){
			if (verbose) print(paste("can't calculate new coeff, (i=",i,") is empty, keep the old one"))
			nw <- w
			nb <- b
			break
		} else{ 
			nw = nwb[-1]
			nb = nwb[1]
		}
		diff = sqrt( sum( (nwb - c(b,w))^2 ) )
		ind = abs(nw) > 0.001  # oracle threshold
		
		#print(table(ind))
		#print(diff)
		
		# if the new model is non-empty, update the old model 
		if (sum(!is.na(ind))>0  & sum(ind)>0) {
			# update w and b, goto the next iteration
			w = nw[ind]
			names(w)=colnames(xtrain)[ind]
			xtrain = xtrain[, ind, drop=FALSE]
			xind = xind[ind]
			b = nb
		}else {
			if (verbose) print(paste("the model in the next iteration, (i=",i,") is empty, keep the old one"))
			diff=tol/2
		}
				
		##  debugging
		#hist(w, main=i)
		#legend("topleft", c(paste("features:", length(w)), 
		#										paste("diff:", round(diff, 5)), 
		#										paste("tol:", tol),
		#										paste("diff - tol: ", round(diff - tol,5) ))  ) 
		##  debugging
		
		i<-i+1
	} # end of while
	if (verbose)print(paste("Elastic SCAD converged in", i-1, "iterations" )) 
	#
		
	ind = abs(nw) > 0.001 # oracle threshold
		
	if (sum(!is.na(ind))>0  & sum(ind)>0){
				
		f = as.vector(xtrain %*% w + b)
		# reduce the calculated time and matrix size
		qx<- .calc.mult.inv_Q_mat2(U=t(x1),D_vec=D_vec, A_vec=Q1_vec, mat2=t(x1) )
		xqx =  0.5 * x1 %*% qx
			
		# Output:
		#   w : direction vector of separating hyperplane
		#   b : the 'bias'
		#   xind : Indices of remained variables
		#   fitted : Fit of regularized SVM (for all patients with reduced set of genes )
		ret <- list(w=w, b=b, xind=xind, index=index, xqx=xqx, fitted=f, type=type, lambda1 = lambda1, lambda2=lambda2, iter=i-1, maxIter=maxIter)
		
		class(ret) <- "scadsvm"
			return(ret)
	} else {
	return("No variable selected.")
	}
}



.calc.test.error.scad<-function(test.x,test.y, f, sg=NULL){
	# calculate test error, confusion table, sensitivity and specificity for test data
	# input:
	## test.x: matrx of ratios, rows - samples, columns - clones
	# test.y - vector of classes for samples
	# f -  model with w,gamma, xind
	# sg =  names of significant (informative) clones
	
	# if we have a model, calulate the test error...
	if (!is.null(f)){
	
		# sepatating line: x'w + b = 0
		sep = as.matrix(test.x)[,f$xind, drop=FALSE] %*% f$w + f$b
		# for classes -1 and 1 :  class = 2*(sep > 0) - 1
		
		class = 2* as.numeric (sep > 0) -1
		
		# 2 different data possible: for GENES and for CLASSES (here only for classes)
				
		# missclassification table for CLASSES		
		tab.classes<-table(class, test.y)
		
		# 3. sensitivity and specificity  for CLASSES
		if (!is.null(tab.classes)){
			if (!nrow(tab.classes)== ncol(tab.classes)) 	  tab.classes<-.extend.to.quad.matrix (tab.classes)
			# sensitivity = TP/ all P = TP /(TP + FN)
			sensitivity.classes<- tab.classes[2,2]/sum(tab.classes[,2])
			# specificity = TN/ all N = TN /(TN + FP)
			specificity.classes <- tab.classes[1,1]/sum(tab.classes[,1])
			# sum von nebendiagonal
	  	# secondary diagonal
			sec.diag<-c(); 	for (j in 1:ncol( tab.classes)) sec.diag<-  c(sec.diag,  tab.classes[ncol( tab.classes)-j+1,j]  )
			error.classes<- ( sum(sec.diag) ) / sum( tab.classes)
		} else {
				error.classes<- NA
				sensitivity.classes<- NA
				specificity.classes<- NA
		}	
	} else{# end of non-empty model f	
		tab.classes=NULL
		error.classes<- NA
		sensitivity.classes<- NA
		specificity.classes<- NA
	}
	
	return(list(classes=list(
								pred.class=class,
								tab=tab.classes,
								error=error.classes,
								sensitivity=sensitivity.classes,
								specificity=specificity.classes )))
	

} 


		