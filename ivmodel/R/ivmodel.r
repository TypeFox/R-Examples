### ivmodel: Generates an instance of ivmodel
###          In this package, our IV model is
###          Y_i = beta_0 + D_i * beta_1 + X_i^T * gamma + epsilon_i
###          E(epsilon_i | Z_i, X_i) = 0
### INPUT: Y, outcome (n * 1 vector)
###        D, exposure (n * 1 vector)
###        Z, instruments (n * L matrix)
###        X, baseline and exogenous covariates (n * p matrix)
###        intercept, should we include the intercept term?
###        beta0: null value for beta
###        alpha: 1 - alpha confidence interval
###        k: values for the k-class estimator
###        heteroSE: heteroskedastic SE?
###        clusterSE: clustered SE?
###        clusterID, if clusterSE is TRUE, then one must provide a vector of
###                   length that's identical to the sample size.
###                   For example, if n = 6 and clusterID = c(1,1,1,2,2,2),
###                   there would be two clusters where the first cluster
###                   is formed by the first three observations and the
###                   second cluster is formed by the last three observations
###                   clusterID can be numeric, character, or factor.
###        deltarange: for sensitivity analysis.
### OUTPUT: ivmodel object which contains a clean-up version of Y,D,Z,X (if available)
ivmodel <- function(Y,D,Z,X,intercept=TRUE,
                    beta0=0,alpha=0.05,k=c(0,1), 
                    heteroSE = FALSE, clusterID = NULL, 
                    deltarange=NULL) {
    
  # Error checking: check to see if necessary inputs are there!
  if(missing(Y)) stop("Y is missing!")
  if(missing(D)) stop("D is missing!")
  if(missing(Z)) stop("Z is missing!")
  
  # Error checking: check Y and D
  if( (!is.vector(Y) && !is.matrix(Y) && !is.data.frame(Y)) || (is.matrix(Y) && ncol(Y) != 1) || (is.data.frame(Y) && ncol(Y) != 1) || (!is.numeric(Y))) stop("Y is not a numeric vector.")
  if( (!is.vector(D) && !is.matrix(D) && !is.data.frame(Y)) || (is.matrix(D) && ncol(D) != 1) || (is.data.frame(D) && ncol(D) != 1) || (!is.numeric(D))) stop("D is not a numeric vector.")
  Y = as.numeric(Y); D = as.numeric(D)
  if(length(Y) != length(D)) stop("Dimension of Y and D are not the same!")
  
  # Error checking: check Z and convert "strings" into factors
  Z = data.frame(Z); stringIndex = sapply(Z,is.character); Z[stringIndex] = lapply(Z[stringIndex],as.factor)
  if(nrow(Z) != length(Y)) stop("Row dimension of Z and Y are not equal!")
  colnames(Z) = paste("Z",colnames(Z),sep="")
  
  # Add intercept as X
  if(intercept && !missing(X)) X = data.frame(X,1)
  if(intercept && missing(X)) X = data.frame(rep(1,length(Y)))
  
  # Error checking: check X and convert "strings" into factors
  if(!missing(X)) {
    X = data.frame(X); stringIndex = sapply(X,is.character); X[stringIndex] = lapply(X[stringIndex],as.factor)
	if(nrow(X) != length(Y)) stop("Row dimension of X and Y are not equal!")
	colnames(X) = paste("X",colnames(X),sep="")	
  }
  
  # Coalesce all data into one data.frame 
  if(!missing(X)) {
    allDataOrig = cbind(Y,D,Z,X)
  } else {
    allDataOrig = cbind(Y,D,Z)
  }
  
  # Fit adjustment model
  ff = terms(Y ~ D + . -1,data=allDataOrig)
  mf = model.frame(ff,allDataOrig)
  
  # Declare Y and D
  Y = as.matrix(mf$Y); colnames(Y) = "Y"
  D = as.matrix(mf$D); colnames(D) = "D"
 
  # Extract Z and X
  allData = sparse.model.matrix(ff,allDataOrig); attr(allData,"assign") = NULL; attr(allData,"contrasts") = NULL
  Zindex = grep(paste("^",colnames(Z),sep="",collapse="|"),colnames(allData))
  Z = allData[,Zindex,drop=FALSE]; colnames(Z) = sub("^Z","",colnames(Z))
  if(!missing(X)) {
    Xindex = grep(paste("^",colnames(X),sep="",collapse="|"),colnames(allData)); X = allData[,Xindex,drop=FALSE]; colnames(X) = sub("^X","",colnames(X))
  }
  # Check to see if there's enough data points after removing missing values
  if(nrow(Y) <= 1) stop("Too much missing data!")

  n = length(Y)
  if(!missing(X)){
    ### clean X and project Y, D, Z
    qrX<-qr(X)
    p = qrrank(qrX)
    if(p==0)
      stop("vector in X are all 0")
    #if(p<ncol(X))
      #X<-qr.Q(qrX)[, 1:p]%*%qrRM(qrX)[1:p, 1:p]
    Yadj = as.matrix(qr.resid(qrX,Y)); Dadj = as.matrix(qr.resid(qrX,D)); Zadj = as.matrix(qr.resid(qrX,Z))
    ### clean Zadj
    ZadjQR = qr(Zadj)
    L = qrrank(ZadjQR)
    if(L==0)
      stop("No useful instrumental variables")
    if(L<ncol(Z)){
      Zadj<-qr.Q(ZadjQR)[, 1:L]%*%qrRM(ZadjQR)[1:L, 1:L]  ### shall we update the Z by doing qr(X, Z)?   
      #qrXZ<-qr(cbind(X, Z))
      #Z<-(qr.Q(qrXZ)[, 1:(p+L)]%*%qrRM(qrXZ)[1:(p+L), 1:(p+L)])[,(p+1):(p+L)]
    }
    ivmodelObject = list(call = match.call(),n=n,L=L,p=p,Y=Y,D=D,Z=Z,X=X,Yadj=Yadj,Dadj=Dadj,Zadj=Zadj, ZadjQR = ZadjQR)

  }else{
    p = 0
    Yadj = Y; Dadj = D; Zadj = as.matrix(Z)
    ### only need to clean Z
    ZadjQR = qr(Zadj)
    L = qrrank(ZadjQR)
    if(L==0)
      stop("No useful instrumental variables")
    if(L<ncol(Z))
      Z<-Zadj<-qr.Q(ZadjQR)[, 1:L]%*%qrRM(ZadjQR)[1:L, 1:L]
	
    ivmodelObject = list(call = match.call(),n=n,L=L,p=p,Y=Y,D=D,Z=Z,X=NA,Yadj=Yadj,Dadj=Dadj,Zadj=Zadj, ZadjQR = ZadjQR)
  }

  class(ivmodelObject) = "ivmodel"

  ivmodelObject$alpha = alpha
  ivmodelObject$beta0 = beta0
  
  ivmodelObject$AR = AR.test(ivmodelObject,beta0=beta0,alpha=alpha)
  ivmodelObject$ARsens = ARsens.test(ivmodelObject,beta0=beta0,alpha=alpha,deltarange=deltarange)
  ivmodelObject$kClass = KClass(ivmodelObject,beta0=beta0,alpha=alpha,k=k,heteroSE=heteroSE,clusterID=clusterID)
  ivmodelObject$LIML = LIML(ivmodelObject,beta0=beta0,alpha=alpha,heteroSE=heteroSE,clusterID=clusterID)  
  ivmodelObject$Fuller = Fuller(ivmodelObject,beta0=beta0,alpha=alpha,heteroSE=heteroSE,clusterID=clusterID)  
  ivmodelObject$CLR = CLR(ivmodelObject,beta0=beta0,alpha=alpha)

  return(ivmodelObject)
}