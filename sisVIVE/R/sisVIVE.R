sisVIVE <- function(Y,D,Z,intercept=TRUE,normalize=TRUE) {
  # Error checking
  if( (!is.vector(Y) && !is.matrix(Y)) | !is.numeric(Y)) stop("Y is not a numeric vector.")
  if( (!is.vector(D) && !is.matrix(D)) | !is.numeric(D)) stop("D is not a numeric vector.")
  if( !is.matrix(Z) | !is.numeric(Z)) stop("Z is not a numerical matrix.")
  if(nrow(Z) != length(Y)) stop("The dimension of Y differs from the row dimension of Z")
  if(nrow(Z) != length(D)) stop("The dimension of D differs from the row dimension of Z")
  if(length(D) != length(Y)) stop("The dimension of Y differs from the dimension of D")
  if(intercept) {
    if( (ncol(Z) + 1) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  } else {
    if( ncol(Z) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  }
  
  # Define constants 
  n = nrow(Z); L = ncol(Z);
  if(intercept) {
    meanY = mean(Y); meanD = mean(D); meanZ = colMeans(Z);
	Y = Y - mean(Y); D = D - meanD; Z = scale(Z,center=TRUE,scale=FALSE);
  } 
  if(normalize) {
  	normZ = sqrt(colSums(Z^2))
    Z = scale(Z,center=FALSE,scale=TRUE) / sqrt(n-1) 
  } else {
  	normZ = rep(1,L)
  }
  
  QR = qr(Z); Yhat = qr.fitted(QR,Y); Dhat = qr.fitted(QR,D);
  Z.transformed = Z - Dhat %*% t(Dhat) %*% Z / sum( Dhat^2)
  Y.transformed = Yhat - Dhat * (sum(Dhat * Yhat) / sum( Dhat^2))
  fit = lars(Z.transformed,Y.transformed,intercept = FALSE,normalize=TRUE)
  
  alpha = coef(fit) # this is unstandardized from lars, but with respect to potentially standardized Z
  alphaSuppSize = apply(alpha,1,function(x){sum(x != 0)}); indexProperSupp = which(alphaSuppSize < L); 
  beta = drop(t(Dhat) %*% (as.numeric(Y) - Z %*% t(alpha) )) / sum(Dhat^2)
  alpha = scale(alpha,FALSE,normZ) # brings alpha back to original Z scale
  attr(alpha,"scaled:scale") = NULL
  
  # Packaging Object
  lambda = fit$lambda;
  whichInvalid = apply(alpha,1,function(x){paste(which(abs(x) > 0),sep="",collapse=",")})
  dimnames(alpha) = list(1:length(beta),dimnames(Z)[[2]]);
  names(beta) = paste(1:length(beta));
  names(whichInvalid) = paste(1:length(beta))
  
  object <- list(call = match.call(), alpha = alpha, beta = beta, whichInvalid = whichInvalid,nInvalid = alphaSuppSize,
                 lambda = lambda,larsObject = fit,Y=Y,D=D,Z=Z,Dhat=Dhat,normZ = normZ)
  class(object) <- "sisVIVE"
  return(object)
}

print.sisVIVE <- function(x,...) {
  object = x
  cat("\nCall:\n")
  dput(object$call)
  printOut = data.frame(object$beta,object$nInvalid, object$whichInvalid)
  colnames(printOut) = c("Estimates of Beta",
                          "   Number of Invalid IVs","Invalid IVs")
  print(printOut,row.names=FALSE)
  invisible(object)
}

predict.sisVIVE <- function(object,lambda,type=c("coefficients","instruments"),...) {
  type = match.arg(type)
  if(!(type %in% c("coefficients","instruments"))) stop("For type, specify either `coefficients' or `instruments'")
  if(missing(lambda)) lambda = object$lambda
  alpha = predict(object$larsObject,s = lambda,type="coefficients",mode="lambda",...)$coefficient
  if(is.vector(alpha)) alpha = matrix(alpha,1,length(alpha))
  if(type == "coefficients") {
    beta = drop(t(object$Dhat) %*% (as.numeric(object$Y) - object$Z %*% t(alpha)))/sum(object$Dhat^2) 
    alpha = scale(alpha,FALSE,object$normZ) # brings alpha back to original Z scale
    attr(alpha,"scaled:scale") = NULL
	return(list(lambda = lambda,alpha = alpha,beta=beta))
  } 
  if(type == "instruments") {
    return(list(lambda = lambda, instruments = apply(alpha,1,function(x){paste(which(abs(x) > 0),sep="",collapse=",")})))
  } 		
}

summary.sisVIVE <- function(object,...) {
  print.sisVIVE(object,...)
}

cv.sisVIVE <- function(Y,D,Z,lambdaSeq,K = 10,intercept=TRUE,normalize=TRUE) {
  # Error checking
  if( (!is.vector(Y) && !is.matrix(Y)) | !is.numeric(Y)) stop("Y is not a numeric vector.")
  if( (!is.vector(D) && !is.matrix(D)) | !is.numeric(D)) stop("D is not a numeric vector.")
  if( !is.matrix(Z) | !is.numeric(Z)) stop("Z is not a numerical matrix.")
  if(nrow(Z) != length(Y)) stop("The dimension of Y differs from the row dimension of Z")
  if(nrow(Z) != length(D)) stop("The dimension of D differs from the row dimension of Z")
  if(length(D) != length(Y)) stop("The dimension of Y differs from the dimension of D")
  if(intercept) {
    if( (ncol(Z) + 1) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  } else {
    if( ncol(Z) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  }
  if(!is.numeric(K) | K > length(Y)) stop("K is not a proper numeric number")
  fitall = sisVIVE(Y,D,Z,intercept=intercept,normalize=normalize)
  if(missing(lambdaSeq) || all(is.na(lambdaSeq))) {
    #warning("Lambda sequence not provided; defaulting to using lambdas provided by sisVIVE")
    lambdaSeq = c(fitall$lambda,seq(from=min(fitall$lambda,na.rm=TRUE),to=2*max(fitall$lambda,na.rm=TRUE),length.out = 100))
    lambdaSeq = sort(unique(lambdaSeq))
  }
  if(any(is.na(lambdaSeq))) {
    warning("Some lambda values are missing. Ignoring these lambda values for cross-validation")
	lambdaSeq = lambdaSeq[!is.na(lambdaSeq)]
  }
  if(length(lambdaSeq) < 2) stop("Only one lambda provided. Please provide multiple lambdas")
  lambdaSeq = sort(lambdaSeq,decreasing=TRUE)
  
  # Define constants
  n = nrow(Z); L = ncol(Z);
 
  # Cross validation
  sampleIndex = sample(rep(1:K,length.out= n))
  errormat = matrix(0,K,length(lambdaSeq))
  for(i in 1:K) {
    testSet = (sampleIndex == i)
	trainfit = sisVIVE(Y[-testSet], D[-testSet], Z[-testSet, , drop=FALSE], intercept = intercept,normalize = normalize)
	trainfitCoef = predict(trainfit, lambda = lambdaSeq, type = "coefficients")
	Y.test = Y[testSet]; D.test = D[testSet]; Z.test = Z[testSet,,drop=FALSE]
	if(intercept) {
      meanY.test = mean(Y.test); meanD.test = mean(D.test); meanZ.test = colMeans(Z.test)
	  Y.test = Y.test - meanY.test; D.test = D.test - meanD.test; Z.test = scale(Z.test,center=TRUE,scale=FALSE)
    } 
	
	QR = qr(Z.test)
	residTest = (as.numeric(Y.test) - Z.test %*% t(trainfitCoef$alpha) - D.test %*% t(trainfitCoef$beta))
	errormat[i,] = colSums(qr.fitted(QR, residTest)^2) #does (PZ %*% residual)^2 summed across observations
  }
  cv = colMeans(errormat)
   
  if(all(is.nan(cv))) {
    warning("All lambdas were invalid. Please try different values of lambda for cross-validation to work")
	return(list(lambda = rep(NA,length(lambdaSeq)),estCVError = NA, alpha = rep(NA,ncol(Z)),beta = NA, whichInvalid = NA))
  } else {
    stderror = apply(errormat,2,function(x){sd(x)/sqrt(K)})
    mincv.index = which.min(cv); onestderrorbound = cv[mincv.index] + stderror[mincv.index]
    onestdLambda = max(lambdaSeq[which( (cv <= onestderrorbound) & (cv >= cv[mincv.index]))]) #this is never empty vector 
    returnOut = predict(fitall,onestdLambda,type="coefficients")
    return(list(lambda = returnOut$lambdaSeq, estCVError = cv[which(onestdLambda == lambdaSeq)], alpha = drop(returnOut$alpha), beta = drop(returnOut$beta),whichInvalid = paste(which(abs(returnOut$alpha) > 0),sep="",collapse=",")))
  }
}
