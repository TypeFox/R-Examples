#library("mixture")
#library("mnormt")

############################################
## Classification matrix given the groups ##
############################################

class <- function(groups,k){
  
  n <- length(groups)
  z <- array(0,c(n,k),dimnames=list(1:n,paste("comp",1:k,sep="")))
  for(i in 1:n) 
    z[i,groups[i]] <- 1
  return(z)
  
}

#################################################################
## Modified Weighted Covariance Matrix for contaminated models ##
#################################################################

.cov.wt.punzo <- function(x,wt,fact){
  
  # fact contains the corrections factors due to the contamination 
  
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  else if (!is.matrix(x)) 
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x))) 
    stop("'x' must contain finite values only")
  p   <- ncol(x)
  n   <- nrow(x)
  mu  <- array(colSums(wt*fact/sum(wt*fact)*x),c(p),dimnames=list(paste("X.",1:p,sep="")))
  cov <- array(crossprod(sqrt(wt*fact/sum(wt))*(x-matrix(rep(mu,n),n,p,byrow=TRUE))),c(p,p),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep="")))
  return(
    list(
      center = mu,
      cov = cov      
    )
  )
}

################################
## M-step of the EM algorithm ##
################################

.m.step <- function(data=NULL, covtype=NULL, w=NULL, fact=matrix(1,nrow(w),ncol(w)), v=1, D=NULL, mtol=NULL, mmax=NULL) {
  G= ncol(w);
  d= ncol(data);
  Sk = array(0, c(d,d,G) )
  gpar= list()
  for (k in 1:G ) {
    gpar[[k]] = list()  	
    temp = .cov.wt.punzo(x=data, wt=w[,k], fact=fact[,k])
    gpar[[k]]$mu    = temp$center
    if (!any(is.na( temp$cov))) gpar[[k]]$sigma = temp$cov
    Sk[,,k]   = temp$cov
  }
  gpar$pi = apply(w,2,mean)
  
  temp = .model.type(modelname = covtype, Sk=Sk, ng=gpar$pi, D=D, mtol= mtol, mmax= mmax )
  gpar$D = temp$D
  for (k in 1:G ) {
    gpar[[k]]$sigma    = temp$sigma[,,k] #+ diag( 1-v , d, d)
    gpar[[k]]$invSigma = temp$invSigma[,,k] 
    gpar[[k]]$logdet   = temp$logdet[k]
  }
  return(gpar)
}

.npar.model <- function(modelname=NULL, p=NULL, G=NULL) {
  val = numeric(3)
  val[1] = G-1 
  val[2] = G*p
  val[3] = .ncovpar(modelname= modelname, p=p, G=G)
  val = sum(val)
  return(val)
}

.ncovpar <- function(modelname=NULL, p=NULL, G=NULL) {
  if (is.null(p)) stop("p is null")
  if (is.null(G)) stop("G is null")
  if (is.null(modelname)) stop("modelname is null")
  
       if (modelname == "EII") npar = 1
  else if (modelname == "VII") npar = G
  else if (modelname == "EEI") npar = p
  else if (modelname == "VEI") npar = p + G -1	
  else if (modelname == "EVI") npar = p*G - G +1
  else if (modelname == "VVI") npar = p*G
  else if (modelname == "EEE") npar = p*(p+1)/2
  else if (modelname == "EEV") npar = G*p*(p+1)/2 - (G-1)*p	
  else if (modelname == "VEV") npar = G*p*(p+1)/2 - (G-1)*(p-1)
  else if (modelname == "VVV") npar = G*p*(p+1)/2
  else if (modelname == "EVE") npar = p*(p+1)/2 + (G-1)*(p-1)
  else if (modelname == "VVE") npar = p*(p+1)/2 + (G-1)*p
  else if (modelname == "VEE") npar = p*(p+1)/2 + (G-1)
  else if (modelname == "EVV") npar = G*p*(p+1)/2 - (G-1)
  else stop("modelname is not correctly defined")
  
  return(npar)		
}

.model.type <- function(modelname=NULL, Sk=NULL, ng=NULL, D=NULL, mtol=1e-10, mmax=10) {
  if (is.null(modelname)) stop("modelname is null")
  
       if (modelname == "EII") val = .msEII(Sk=Sk, ng=ng)
  else if (modelname == "VII") val = .msVII(Sk=Sk, ng=ng)
  else if (modelname == "EEI") val = .msEEI(Sk=Sk, ng=ng)
  else if (modelname == "VEI") val = .msVEI(Sk=Sk, ng=ng, eplison= mtol, max.iter= mmax)
  else if (modelname == "EVI") val = .msEVI(Sk=Sk, ng=ng)
  else if (modelname == "VVI") val = .msVVI(Sk=Sk, ng=ng)
  else if (modelname == "EEE") val = .msEEE(Sk=Sk, ng=ng)
  else if (modelname == "EEV") val = .msEEV(Sk=Sk, ng=ng)
  else if (modelname == "VEV") val = .msVEV(Sk=Sk, ng=ng, eplison= mtol, max.iter= mmax)
  else if (modelname == "VVV") val = .msVVV(Sk=Sk, ng=ng)
  else if (modelname == "EVE") val = .msEVE(Sk=Sk, ng=ng, D0=D, eplison= mtol, max.iter= mmax)
  else if (modelname == "VVE") val = .msVVE(Sk=Sk, ng=ng, D0=D, eplison= mtol, max.iter= mmax)
  else if (modelname == "VEE") val = .msVEE(Sk=Sk, ng=ng, eplison= mtol, max.iter= mmax)
  else if (modelname == "EVV") val = .msEVV(Sk=Sk, ng=ng)
  else stop("modelname or covtype is not correctly defined")
  
  if (!is.list(val)) val = list(sigma=val)
  return(val)		
}

.gpcm <- function(data=NULL,  G=1:3, mnames=NULL, start=0, label=NULL, veo=FALSE, nmax=1000, atol=1e-8, mtol=1e-8, mmax=10, pprogress=FALSE, pwarning=TRUE) {
  if (is.null(data)) stop('Hey, we need some data, please! data is null')
  if (!is.matrix(data)) stop('The data needs to be in matrix form')
  if (!is.numeric(data)) stop('The data is required to be numeric')
  if (nrow(data) == 1) stop('nrow(data) is equal to 1')
  if (ncol(data) == 1) stop('ncol(data) is equal to 1; This function currently only works with multivariate data p > 1')
  if (any(is.na(data))) stop('No NAs allowed.')
  
  if (is.null(G)) stop('G is NULL')
  G = as.integer(ceiling(G))
  if (!is.integer(G)) stop('G is not a integer')
  if ( any(G < 1)) stop('G is not a positive integer')
  
  if (is.null(mnames) )  mnames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV", "EVE", "VVE", "VEE", "EVV")
  
  bic = array(0, dim= c(length(G), length(mnames), 3), dimnames=list(G, mnames, c('loglik', "npar", "BIC")) )
  model = NULL; curBIC = Inf;
  for (g in 1:length(G)) {
    for (i in 1:length(mnames)) {
      if ( pprogress ) print(c(G[g],mnames[i]))
      if (veo | .npar.model(modelname=mnames[i], p=ncol(data), G=G[g]) < nrow(data)) {
        a =	try( { .EM(data=data, G=G[g], nmax=nmax, covtype=mnames[i], start=start, label=label, atol= atol, mtol=mtol, mmax=mmax ) }, pwarning)
        
        if ( length(a) > 1) {
          bic[g,i,1:2] = c(a$loglik[length(a$loglik)], .npar.model(modelname=mnames[i], p=ncol(data), G=G[g]) )
          bic[g,i,3] =  -2*bic[g,i,1] + bic[g,i,2]*log(nrow(data)) 
          if ( bic[g,i,3] < curBIC) {
            model  = append(a,list(mtype=mnames[i]))
            curBIC = bic[g,i,3] 
          }
        } else {
          bic[g,i,1:2] = c(NA, .npar.model(modelname=mnames[i], p=ncol(data), G=G[g]) )
          bic[g,i,3] =  NA
        }
      } else {
        bic[g,i,1:2] = c(NA, .npar.model(modelname=mnames[i], p=ncol(data), G=G[g]) )
        bic[g,i,3] =  NA	
      } 	
      
    }}
  
  if ( is.null(start) ) startobject= "deterministic annealing with default values."
  else if ( is.matrix(start)) startobject= "a user specified initialization  matrix."
  else if ( is.function(start)) startobject= "a user specified initialization function."
  else if ( length(start) > 1 ) startobject= "deterministic annealing with user specified values"
  else if ( start == 0) startobject= "k-means"
  else if ( start > 0) startobject= paste( start, " random initializations", collapse="" )
  
  bicModel = list(G=length(model$gpar$pi), covtype=model$mtype, bic=curBIC )
  val = list( start=start, startobject= startobject, gpar=model$gpar, loglik=model$loglik, z=model$z, map=model$MAP, BIC=bic, bicModel= bicModel)
  #	val = list( model=model, BIC=bic, bicModel= bicModel)
  
  class(val)<-"gpcm"
  return(val)
}

.print.gpcm <-function(x, ...){
  #    a = object$BIC[,,1]
  cat("The model choosen by applying the BIC criteria has ", x$bicModel$G, "component(s) and the ", x$bicModel$covtype, "covariance structure\n using ", x$startobject, "\n"  )
  #endPrint(a)
}
.summary.gpcm <- function(object, ...){
  bicl = object$BIC[,,3]
  cat("BIC for each model, number of components (rows), and covariance structure (columns).\n")
  print.default( bicl)
}

.plot.gpcm <- function(x, ...) {
  bicl = x$BIC[,,1]
  g   = dimnames(bicl)[[1]]
  cov = dimnames(bicl)[[2]]
  ncov = length(cov)
  ng   = length(g)
  
  plot( g, bicl[,1],  ylim=range(bicl,na.rm=TRUE), xlab="# components", ylab="BIC",type='l')
  if (ncov >1) {
    for (i in 2:ncov) lines( g, bicl[,i], col=i)
  }
  legend( "bottomright", legend=cov, col=1:ncov, lty=1)
  
}

.igpar = function(data=NULL, g=NULL, covtype=NULL, start=NULL, labels=NULL, mtol=NULL, mmax=NULL) {
  if (is.null(start)) start = seq( 0.1, 1, length.out=5)
  
  if (is.function(start) ) w = start(data=data, g=g, covtype=covtype)	
  else if ( !is.null(dim(start)) ) w = start
  else if (length(start) > 1) w = .igparv(data=data, g=g, covtype=covtype, vseq0=start , mtol=mtol, mmax = mmax, labels=labels)
  else if (start == 0) w = .igpark(data=data, g=g, covtype=covtype)	
  else if (start >  0) w = .rgpar(data=data, g=g, covtype=covtype, n=ceiling(start), labels=labels, mtol=mtol, mmax = mmax )	
  else stop(paste('Initialization method ', start, " does not compute!"))
  
  if (!is.matrix(w)) stop("The zij initialization matrix is not a matrix")
  if (any(w < 0)) stop("Some of the elements of the zij initialization matrix are less than zero")
  #	if ( any(apply(w,1,sum) !=1)  ) stop("Some of the rows of the zij initialization matrix do not sume to 1")
  if (nrow(w) < nrow(data) ) stop("The nrow(zij) of the initialization matrix is less than nrow(data)")
  if (nrow(w) > nrow(data) ) stop("The nrow(zij) of the initialization matrix is greater than nrow(data)")
  if (ncol(w) < g ) stop("The nrow(zij) of the initialization matrix is less than  match g")
  if (ncol(w) > g ) stop("The nrow(zij) of the initialization matrix is greater than  match g")
  
  w = .combinewk(w, label=labels)
  gpar = .m.step(data=data, covtype=covtype, w=w, v=1, mtol=mtol, mmax = mmax)
  return(gpar)
}

.igparv <- function(data=NULL, g=NULL,covtype=NULL, vseq0=NULL, labels=NULL, mtol=NULL, mmax = NULL) {
  vseq0 = as.numeric(vseq0)
  if (is.null(vseq0)) stop('The sequence for deterministic annealing is NULL')
  if ( !all( vseq0 <=1 &  vseq0 >=0)  ) stop('The sequence for deterministic annealing must be between 0 and 1')
  
  w = .rwgpar(data = data, g=g, covtype = covtype)
  w = .combinewk(w, label=labels)
  gpar = .m.step(data=data, covtype=covtype, w=w, v=1, mtol= mtol, mmax = mmax)
  gpar = .EMv(data=data, gpar0 = gpar, G=g, vseq=vseq0, m=1, label=labels, covtype = covtype, mtol= mtol, mmax = mmax )$gpar
  w    = .weights(data=data, gpar= gpar) 
  return(w)
}

.rgpar <- function(data=NULL, g=NULL, covtype=NULL, n=1, labels=NULL, mtol=NULL, mmax = NULL) {
  w.old      = .rwgpar(data= data, g=g,covtype= covtype, labels=labels)
  gpar.old   = .m.step(data=data, covtype=covtype, w=w.old, v=1, mtol=mtol, mmax = mmax)
  loglik.old = .loglik(data=data, gpar=gpar.old)
  
  for (i in 1:n) { 
    w.new      = .rwgpar(data= data, g=g,covtype= covtype, labels=labels)
    gpar.new   = .m.step(data=data, covtype=covtype, w=w.new, v=1, mtol=mtol, mmax = mmax)
    loglik.new = .loglik(data, gpar.new)
    
    if (loglik.new > loglik.old) {
      w.old      = w.old
      gpar.old   = gpar.new
      loglik.old = loglik.new
    }
  }		
  return(w.old)
}

.rwgpar <- function(data=NULL, g=NULL,covtype=NULL, labels=NULL) {
  w = matrix(rexp(nrow(data)*g),nrow=nrow(data),ncol=g)
  w = matrix(t(apply(w,1, function(z) { z/sum(z)})),nrow=nrow(data),ncol=g)
  w = .combinewk(w, label=labels)
  return(w)
}

.igpark <- function(data=NULL, g=NULL,covtype=NULL) {
  lw = kmeans(data, centers=g, iter.max=10)$cluster
  w  = .combinewk(matrix(0,nrow=nrow(data),ncol=g), label=lw)
  return(w)
}





.EM <- function(data=NULL, gpar0=NULL, G=2, start=1, label=NULL, covtype=NULL, nmax=1000, atol=1e-8, mtol=1e-8, mmax=10 ) {
  val        = list()
  if (is.null(gpar0)) val$gpar = .igpar(data=data, g=G, covtype=covtype, start=start, labels=label, mtol=mtol, mmax=mmax)
  else val$gpar = gpar0	
  val$loglik = numeric(nmax)
  
  val$loglik[1] = .loglik(data=data, gpar=val$gpar)
  tempw         = .e.step(data=data, gpar=val$gpar, labels=label)
  val$gpar      = .m.step(data=data, covtype=covtype, w=tempw, mtol=mtol, mmax=mmax)
  val$loglik[2] = .loglik(data=data, gpar=val$gpar)
  tempw         = .e.step(data=data, gpar=val$gpar, labels=label)
  val$gpar      = .m.step(data=data, covtype=covtype, w=tempw, mtol=mtol, mmax=mmax)
  val$loglik[3] = .loglik(data=data, gpar=val$gpar)
  i  =  3
  
  while ( ( .getall(val$loglik[(i-2):i]) > atol) & (i < (nmax) ) )  {
    i = i+1
    tempw         = .e.step(data=data, gpar=val$gpar, labels=label)
    val$gpar      = .m.step(data=data, covtype=covtype, w=tempw, D=val$gpar$D, mtol=mtol, mmax=mmax)
    val$loglik[i] = .loglik(data=data, gpar=val$gpar)
  }
  
  val$loglik = val$loglik[1:i]
  val$z   = .e.step(data=data, gpar=val$gpar, labels=label) 
  val$MAP = .MAP(data=data, gpar=val$gpar, label=label) 
  return(val)
}






.EMn <- function(data=NULL, gpar0 = NULL, G=2, n=10, label =NULL, covtype="VVV", mtol=1e-8, mmax=10 ) {
  val = list()
  if (is.null(gpar0)) val$gpar = .igpar(data=data, g=G, covtype=covtype, mtol=mtol, mmax=mmax)
  else val$gpar = gpar0	
  val$loglik = numeric(n)
  for (i in 1:n) {
    tempw         = .e.step(data=data, gpar=val$gpar, labels=label)
    val$gpar      = .m.step(data=data, covtype=covtype, w=tempw, mtol=mtol, mmax=mmax)
    val$loglik[i] = .loglik(data=data, gpar=val$gpar)
  }
  #	val$label  = label
  return(val)
}


.EMv <- function(data=NULL, gpar0=NULL, G=3, vseq=c(1,1), m=2, label=NULL, covtype="VVV", mtol=NULL, mmax = NULL  ) {
  val = list()
  if (is.null(gpar0)) val$gpar = .rgpar(data=data, g=G, covtype=covtype, n=1 )
  else val$gpar = gpar0	
  
  val$loglik = numeric(length(vseq)*m)
  count = 1
  for (i in 1:length(vseq)) { for (j in 1:m) {
    tempw         = .e.step(data=data, gpar=val$gpar, labels=label, v=vseq[i])
    val$gpar      = .m.step(data=data, covtype=covtype, w=tempw, D=val$gpar$D, mtol= mtol, mmax = mmax)
    val$loglik[i] = .loglik(data, val$gpar)
    count = count + 1
    #cat("iteration = ", count-1, "\t v= ", vseq[i], "\t loglik = ", val$loglik[count-1], "\t pi=", val$gpar$pi, "\n")
  }}
  return(val)
}



.weights <- function(data=NULL, gpar=NULL, v=1) {
  d = ncol(data)
  G = length(gpar$pi)	
  if (G > 1) {
    zlog = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
    for (k in 1:G ) zlog[,k] = -1/2*mahalanobis(x=data, center=gpar[[k]]$mu, cov=gpar[[k]]$invSigma, inverted=TRUE) -1/2*gpar[[k]]$logdet - d/2 *log(2*pi)
    
    w = t(apply( zlog, 1, function(z,wt,v) { 
      x= exp( v*(z + log(wt)) ) 
      x=x/sum(x);
      return(x) }, wt=gpar$pi,v=v ))
    #		w = t(apply(zlog, 1, function(z,wt,v) { 
    #			x=exp(v*(z + log(wt)) );
    #			if (sum(x)  == 0) x= rep(1,length(x))
    #			x =  x/sum(x)
    #			return( x ) 
    #			}, wt=gpar$pi,v=v ))
  } else w = matrix(1,nrow=nrow(data), ncol=G)
  return(w)
}


.e.step <- function(data=NULL, gpar=NULL, labels=NULL, v=1) {
  w = .weights(data=data, gpar=gpar,v=v)
  if (!is.null(labels)) w = .combinewk(weights=w, label= labels)
  return(w)
}

.loglik = function(data, gpar) {
  # output is a G x nrow(data) matrix
  d = ncol(data)
  G = length(gpar$pi)
  zlog = matrix(0, nrow=nrow(data), ncol=G)
  for (k in 1:G) zlog[,k] = -1/2*mahalanobis(x=data, center=gpar[[k]]$mu, cov=gpar[[k]]$invSigma, inverted=TRUE) -1/2*gpar[[k]]$logdet - d/2*(log(2)+log(pi))
  
  w = apply( exp(zlog),1,function(z,wt) { sum(z*wt) } , wt=gpar$pi)
  val = sum(log(w))
  if( is.nan(val) ) val =NA	
  return(val)
}



.getall <- function(loglik) {
  lm1 = loglik[3]
  lm  = loglik[2]
  lm_1  = loglik[1]
  am = (lm1 - lm)/(lm - lm_1)
  lm1.Inf = lm + (lm1 - lm)/(1-am)
  val = lm1.Inf - lm	
  if (is.nan(val) ) val =0
  return( abs(val) )
}



.MAP <- function(data, gpar, label=NULL) {
  w = .weights(data=data, gpar=gpar, v=1)
  if (!is.null(label)) w = .combinewk(weights=w, label= label)
  z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
  z = as.numeric(z)
  return( z)	
}


.combinewk <- function(weights=NULL, label=NULL)	{
  # known is a numeric with 
  # 0 if unknown group membership 
  # 1,2,3,.. for label of known group
  if (!is.null(label)) { #stop('label is null')
    label = as.integer(label)
    
    if (any(!is.integer(label)))	stop("Labels are not integers")
    if (any(label < 0 ))	stop("Labels can only be positive integers")
    if (ncol(weights) < max(label) ) stop("Number of groups is less then the number groups given by labels")
    
    if ( sum(label!=0) == nrow(weights) ) {
      if (ncol(weights) > max(label) ) stop("Every observations has a label; Cannot fit more groups to the data then given the by the labels.")
    }
    
    kw     = label !=0
    for (j in 1:ncol(weights)) weights[kw,j] = (label == j)[kw]
  }
  return(weights)	
}






.msEEE <- function(Sk=NULL, ng=NULL) {
  # Sk is an array of with dim (p x p x G)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  W = .sumSk.wt(Sk=Sk, wt=ng, d=d, G=G)/sum(ng)
  
  val = array(0, c(d,d,G))
  for (g in 1:G) val[,,g] = W
  
  logdetW = log(det(W))
  invW    = solve(W)
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = W
    val$invSigma[,,g] = invW
    val$logdet[g]     = logdetW
  }
  return(val)
}

.sumSk.wt <- function(Sk=NULL, wt=NULL, d=NULL, G=NULL) {
  # Sum Sk over the groups with weights wt.
  W = matrix(0, nrow=d, ncol=d )
  for (g in 1:G) W = W + Sk[,,g]* wt[g]
  return(W)
}


.msEEV <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 100) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  EWk = Sk
  A  = matrix(0,d,d)
  for (g in 1:G) {
    Wk = Sk[,,g]*ng[g]
    EWk[,,g] = eigen(Wk)$vectors
    A = A + t(EWk[,,g]) %*% Wk %*% EWk[,,g]
  }
  
  lam = prod(diag(A))^(1/d)
  A   = A/lam
  lam = lam/sum(ng)
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = lam *( EWk[,,g] %*% A %*% t(EWk[,,g]) )
    val$invSigma[,,g] = 1/lam *( EWk[,,g] %*% diag(1/diag(A),d) %*% t(EWk[,,g]) )
    val$logdet[g]     = d*log(lam) 
  }
  return(val)
}

.getA <- function(Ok=NULL, lam=NULL, G=NULL, d=NULL) {
  A = matrix(0,d,d)
  for (g in 1:G) A = A + Ok[,,g]/lam[g]
  A = diag(A)
  A = diag( A/prod(A)^(1/d) )
  return( A )	
}


.getOk <- function(Sk=NULL, ng=NULL, G=NULL) {
  Ok  = Sk
  for (g in 1:G) {
    Wk  = Sk[,,g]*ng[g]
    EWk = eigen(Wk)$vectors
    Ok[,,g] = t(EWk) %*% Wk %*% EWk
    Ok[,,g] = diag(diag(Ok[,,g]))
  }
  return(Ok)	
}

.getEkOk <- function(Sk=NULL, ng=NULL, G=NULL) {
  Ok = Sk
  EWk = Sk
  for (g in 1:G) {
    Wk       = Sk[,,g]*ng[g]
    EWk[,,g] = eigen(Wk)$vectors
    Ok[,,g]  = t(EWk[,,g]) %*% Wk %*% EWk[,,g]
  }
  return(list(Ok=Ok,EWk=EWk))	
}

.msVEV <- function(Sk=NULL, ng=NULL, eplison=1e-14, max.iter= 100) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  temp = .getEkOk(Sk=Sk, ng=ng, G=G)
  Ok  = temp$Ok
  EWk = temp$EWk
  lam = apply(Ok,3,function(z) { sum(diag(z)) } )/(ng*d)
  A   = .getA(Ok=Ok, lam=lam, d=d, G=G) 
  lam = apply(Ok,3,function(z, invA) { sum(diag(z*invA)) }, invA=diag(1/diag(A)) )/(ng*d)
  
  conv = c( d*sum(ng*(1+log(lam))), Inf  )
  count = 1 
  while ( diff(conv)/conv[1] > eplison & count < max.iter) {
    A   = .getA(Ok=Ok, lam=lam, d=d, G=G) 
    lam = apply(Ok,3,function(z, invA) { sum(diag(z*invA)) }, invA=diag(1/diag(A)) )/(ng*d)
    conv = c(d*sum(ng*(1+log(lam))), conv[1] )
    count = count +1
  }
  
  val = array(0, c(d,d,G))
  for (g in 1:G) val[,,g] = lam[g] * ( EWk[,,g] %*% A %*% t(EWk[,,g]) )
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    =  lam[g] * ( EWk[,,g] %*% A %*% t(EWk[,,g]) )
    val$invSigma[,,g] =  1/lam[g] * ( EWk[,,g] %*% diag(1/diag(A),d) %*% t(EWk[,,g]) )
    val$logdet[g]     =  d*log(lam[g])
  }
  return(val)
}

.msVVV <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    =  Sk[,,g]
    val$invSigma[,,g] =  solve(Sk[,,g])
    val$logdet[g]     =  log(det(Sk[,,g]) )
  }
  return(val)
}

.msEEI <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  W = .sumSk.wt(Sk=Sk, wt=ng, d=d, G=G)/sum(ng)
  B = diag(diag(W))
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G) )
  for (g in 1:G) { 
    val$sigma[,,g]    = B
    val$invSigma[,,g] = diag( 1/diag(B),d )
    val$logdet[g]     = sum(log( diag(B) ))
  }
  return(val)
}

.msVEI <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 100) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  lam = apply(Sk,3,function(z) { sum(diag(z)) } )/d
  #	W = sumSk.wt(Sk=Sk, wt=ng*lam, d=d, G=G)
  W = .sumSk.wt(Sk=Sk, wt=ng/lam, d=d, G=G)
  W = diag(W)
  B = diag(W/prod(W)^(1/d))
  lam = apply(Sk,3,function(z, invB) { sum(diag(z*invB)) }, invB=diag(1/diag(B) ) )/d
  
  conv = c( d*sum(ng*(1+log(lam))), Inf  )
  count = 1 
  while ( abs(diff(conv)) > eplison & count < max.iter) {
    #	while ( count < max.iter) {
    #		W = sumSk.wt(Sk=Sk, wt=ng*lam, d=d, G=G)
    W = .sumSk.wt(Sk=Sk, wt=ng/lam, d=d, G=G)
    
    W = diag(W)
    B = diag(W/prod(W)^(1/d))
    lam = apply(Sk,3,function(z, invB) { sum(diag(z*invB)) }, invB=diag(1/diag(B) ) )/d
    
    conv = c(d*sum(ng*(1+log(lam))), conv[1] )
    count = count +1
  }
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = lam[g]*B
    val$invSigma[,,g] = diag( 1/diag(B) * 1/lam[g],d )
    val$logdet[g]     = d*log(lam[g])
  }
  return(val)
}


.msEVI <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  Bk = matrix(0,nrow=d, ncol=G)
  for (g in 1:G) Bk[,g] = diag(Sk[,,g]*ng[g])
  lam = apply(Bk, 2, prod)^(1/d) 
  Bk  = sweep(Bk, 2, 1/lam, FUN="*")	
  lam = sum(lam)/sum(ng)
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = lam * diag(Bk[,g],d)
    val$invSigma[,,g] = 1/lam * diag(1/Bk[,g],d)
    val$logdet[g]     = d*log(lam)
  }
  return(val)
}




.msVVI <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  #	Bk = matrix(0,nrow=d, ncol=G)
  #	for (g in 1:G) Bk[,g] = diag(Sk[,,g]*ng[g])
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = diag(diag(Sk[,,g]),d)
    val$invSigma[,,g] = diag(1/diag(Sk[,,g]),d)
    val$logdet[g]     = sum( log(diag(Sk[,,g])) )
  }
  return(val)
}

.msEII <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  W = .sumSk.wt(Sk=Sk, wt=ng, d=d, G=G)/sum(ng)
  lam = sum(diag(W))/(sum(ng)*d)
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G) )
  for (g in 1:G) { 
    val$sigma[,,g]    = diag(rep(lam,d), d)
    val$invSigma[,,g] = diag(rep(1/lam,d),d)
    val$logdet[g]     = d*log(lam) 
  }
  return(val)
}


.msVII <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    sumdiagSkg        = sum(diag(Sk[,,g]))
    val$sigma[,,g]    = diag(rep(sumdiagSkg/d,d))
    val$invSigma[,,g] = diag(rep( d/sumdiagSkg, d))
    val$logdet[g]     = d* log( sumdiagSkg ) - d*log(d)
  }
  return(val)
}

.msVVE <- function(Sk=NULL, ng=NULL) {
  # Sk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  #	lam = numeric(G)
  #	for (g in 1:G) 	lam[g] = sum(diag(Sk[,,g]))/d
  
  val = array(0, c(d,d,G))
  for (g in 1:G) val[,,g] = diag(rep(sum(diag(Sk[,,g]))/d,d))
  return(val)
}



##################### 
### new.cov models
##################### 

.newD3.MM <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL, tmax = 100) {
  z = matrix(0,d,d)
  lambda =0 
  for (g in 1:G) {
    lambdak = max(eigen(Wk[,,g])$values)
    z = z + diag(1/Ak[,g]) %*% t(D) %*% Wk[,,g]  -lambdak *(diag(1/Ak[,g])%*% t(D) )
  } 
  z1 = svd(z)
  Xk1 = (z1$v) %*% t(z1$u) 
  return( Xk1 )
}

.newD4.MM <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL, tmax = 100) {
  z = matrix(0,d,d)
  lambda =0 
  for (g in 1:G) {
    lambdak = max(1/Ak[,g])
    #z = z + diag(1/Ak[,g]) %*% t(D) %*% Wk[,,g]  -lambdak *(diag(1/Ak[,g])%*% t(D) )
    z = z + Wk[,,g] %*% (D) %*% diag(1/Ak[,g])  -lambdak *(Wk[,,g]%*% (D) )
  } 
  z1 = svd(z)
  #	Xk1 = (z1$v) %*% t(z1$u) 
  #	return( t(Xk1) )
  # OR 
  Xk1 = (z1$v) %*% t(z1$u) 
  return( t(Xk1) )
  
}

.newD <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL, tmax = 100) {
  D6 = D
  D6 = .newD3.MM(D=D6, d=d, G=G, Wk=Wk, Ak=Ak, tmax = 100)
  D6 = .newD4.MM(D=D6, d=d, G=G, Wk=Wk, Ak=Ak, tmax = 100)
  return(D6)
}

.testval <- function(Wk=NULL, Ak=NULL, D=NULL, G=NULL) {
  z = numeric(G)
  #	for (g in 1:G) z[g] = sum(diag( D %*%  diag(1/Ak[,g]) %*% t(D) %*% Wk[,,g]))
  for (g in 1:G) z[g] = sum(diag( t(D) %*% Wk[,,g] %*% D %*%  diag(1/Ak[,g])  ))
  return(sum(z))
}

.testgrad.D <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL) {
  z = matrix(0,d,d)
  for (g in 1:G) z = z + Wk[,,g] %*% D %*% diag(1/Ak[,g]) 
  return(2*z)
}



.msEVE <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 10, D0=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1]; G = dim(Sk)[3];
  
  Wk = Sk
  W  = matrix(0,d,d)
  for (g in 1:G) {
    Wk[,,g] = Sk[,,g]*ng[g]
    W  = W + Wk[,,g]
  }
  #	D  = diag(rep(1,d))
  D = D0
  if (is.null(D)) D  = diag(rep(1,d))
  #	D  = t(eigen(W)$vectors)
  Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% (D) ) }, D=D )
  Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})
  #print(c(0, 1, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
  D = .newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d)
  #print(c(0, 2, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
  
  conv = c( .testval(Wk=Wk,Ak=Ak,D=D,G=G), Inf  )
  count = 1  
  while ( diff(conv)/abs(conv[1]) > eplison & count < max.iter) {
    D = .newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d) 
    Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% (D) ) }, D=D )
    Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})
    
    #print(c(count, 0, .testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
    conv = c(.testval(Wk=Wk,Ak=Ak,D=D,G=G), conv[1] )
    count = count +1
  }
  lam =  0
  for (g in 1:G) lam = lam  + sum(diag( D %*% diag(1/Ak[,g])%*% t(D) %*% Wk[,,g] ))
  lam = lam/(sum(ng)*d)
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G), D=D  )
  for (g in 1:G) { 
    val$sigma[,,g]    = D %*% diag(lam*Ak[,g])%*% t(D)
    val$invSigma[,,g] = D %*% diag(1/lam*1/Ak[,g])%*% t(D)
    val$logdet[g]     = d*log( lam ) 
  }
  return(val)
  
}


.msVVE <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 10, D0=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1]; G = dim(Sk)[3];
  
  Wk = Sk
  W  = matrix(0,d,d)
  for (g in 1:G) {
    Wk[,,g] = Sk[,,g]*ng[g]
    W  = W + Wk[,,g]
  }
  
  D = D0
  if (is.null(D)) D  = diag(rep(1,d))
  
  Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% D ) }, D=D )
  Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})
  
  #print(c(0, 1, .testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
  D = .newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d)
  #print(c(0, 2, .testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
  
  conv = c( .testval(Wk=Wk,Ak=Ak,D=D,G=G), Inf  )
  count = 1  
  while ( diff(conv)/abs(conv[1]) > eplison & count < max.iter) {
    Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% D ) }, D=D )
    Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})
    D = .newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d) 
    
    conv = c(.testval(Wk=Wk,Ak=Ak,D=D,G=G), conv[1] )
    count = count +1
  }
  #print(count)
  lam = numeric(G) 
  for (g in 1:G) lam[g] =sum(diag( D %*% diag(1/Ak[,g])%*% t(D) %*% Sk[,,g] ))/d
  #print(apply(Ak,2,prod))
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = (D %*% diag(lam[g]*Ak[,g])%*% t(D))
    val$invSigma[,,g] = (D %*% diag(1/lam[g]*1/Ak[,g])%*% t(D))
    val$logdet[g]     =  d*log(lam[g])
  }
  return(val)
}



.msVEE <- function(Sk=NULL, ng=NULL, eplison=1e-14, max.iter=100) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  Wk = Sk
  W  = matrix(0,d,d)
  for (g in 1:G) {
    Wk[,,g] = Sk[,,g]*ng[g]
    W  = W + Wk[,,g]
  }
  
  C = W/det(W)^(1/d)
  invC = solve(C)
  lam = apply(Sk,3,function(z, invC) { sum(diag(z %*% invC)) }, invC=invC)/d
  
  val1 = sum(apply(Wk,3,function(z, invC) { sum(diag(z*invC)) }, invC= invC)/lam) +d*sum(ng*lam)
  conv = c(val1, Inf  )
  count = 1 
  while ( diff(conv)/conv[1]  > eplison & count < max.iter) {
    #for (i in 1:max.iter) {
    C = .sumSk.wt(Sk=Wk, wt=1/lam, d=d,G=G)
    C = C/det(C)^(1/d)		
    invC = solve(C)
    
    lam   = apply(Sk,3,function(z, invC) { sum(diag(z %*% invC)) }, invC=invC)/d
    val1  = sum(apply(Wk,3,function(z, invC) { sum(diag(z*invC)) }, invC= invC )/lam) +d*sum(ng*lam)
    conv  = c(val1, conv[1] )
    count = count +1
  }
  #print(c(count,det(C), lam)	)
  
  invC = solve(C)
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = lam[g]*C
    val$invSigma[,,g] = 1/lam[g]*invC
    val$logdet[g]     = d*log(lam[g])
  }
  return(val)
}


.msEVV <- function(Sk=NULL, ng=NULL, eplison=1e-12, max.iter= 100) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  Wk = Sk
  Ck = Sk
  lam = numeric(G)
  for (g in 1:G) {
    Wk[,,g] = Sk[,,g]*ng[g]
    lam[g]  = det(Wk[,,g])^(1/d)
  }
  Ck = sweep(Wk, 3, 1/lam, FUN="*")	
  lam = sum(lam)
  
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = lam * Ck[,,g]
    val$invSigma[,,g] = 1/lam * solve(Ck[,,g])
    val$logdet[g]     =  d* log(lam ) 
  }
  return(val)
}

####################################
## Mixture of Contaminated Models ##
####################################

.MCM <- function(
	X,			                      # matrix of data
  k,                            # number of groups
	initialization="mclust",      # initialization procedure: "random.soft", "random.hard", "manual", or "mclust" (which is nested)
	modelname="VVV",              # one of the 14 models of Celeaux & Govaert (1995)       
	alphacon=TRUE,                # if TRUE, alpha is constrained to be >0.5
	alphamin=rep(0.5,k),          # minimum proportion of good data in each group 
  alphafix=FALSE,               # if TRUE, the inner weights are not estimated
	alpha=NULL,                   # vector of dimension k with proportion of good observations in each group
  etacon=TRUE,                  # if TRUE, eta is constrained to be >1
	etafix=FALSE,                 # if TRUE, eta is not estimated
	eta=NULL,                     # inflation parameters
  etamax=200,                   # maximum value of eta
	start.z=NULL,                 # (n x k)-matrix of soft or hard classification: it is used only if initialization="manual"		
	start.v=NULL,                 # (n x 2 x k)-array of soft or hard classification in each group: it is used only if initialization="manual"  	
	start=0,                      # initialization for the package mixture
  ind.label=NULL,               # indexes of the labelled observations
  label=NULL,                   # groups of the labelled observations
  iter.max=1000,                # maximum number of iterations in the EM-algorithm
	threshold=1.0e-04             # stopping rule in the Aitken rule
	#loglikplot=TRUE,              # if TRUE, the log-likelihood values against the iterations are plotted
	#plot=TRUE                     # if TRUE, the plot of the classified data is given
	)
{
  
#call = match.call()
  
  if (is.data.frame(X)) 
    X <- as.matrix(X)
    
  n         <- nrow(X)    # sample size
  p         <- ncol(X)    # number of variables
    
  if (is.null(X))     stop('Hey, we need some data, please! X is null')
  if (!is.matrix(X))  stop('X needs to be in matrix form')
  if (!is.numeric(X)) stop('X is required to be numeric')
  if (n == 1)   stop('nrow(X) is equal to 1')
  if (p == 1)   stop('ncol(X) is equal to 1; This function currently only works with multivariate data p > 1')
  if (any(is.na(X)))  stop('No NAs allowed.')
  if (is.null(k)) stop('k is NULL')
  k <- as.integer(ceiling(k))
  if (!is.integer(k)) stop('k is not a integer')
  if (any(k < 1)) stop('k is not a positive integer')

# for model-based classification

lab <- NULL
if(is.vector(label)){
  nlab   <- length(label)
  nunlab <- n-nlab
  lab    <- numeric(n)
  lab[ind.label] <- label
}  

prior     <- numeric(k) # proportion of each group
priorgood <- numeric(k) # proportion of good observations in each group
v         <- array(0.5,c(n,2,k),dimnames=list(1:n,c("good","bad"),paste("group ",1:k,sep="")))
mu        <- array(0,c(p,k),dimnames=list(paste("X.",1:p,sep=""),paste("group ",1:k,sep="")))
Sigma     <- array(0,c(p,p,k),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:k,sep=""))) 
invSigma  <- array(0,c(p,p,k),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:k,sep=""))) 
lambda    <- array(0,c(k),dimnames = list(paste("group ",1:k,sep="")))
Delta     <- array(0,c(p,p,k),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:k,sep="")))
Gamma     <- array(0,c(p,p,k),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:k,sep="")))
if(!etafix)
  eta  <- rep(2,k)  # inflation parameters
correction <- array(0,c(n,k),dimnames=list(1:n,paste("group ",1:k,sep=""))) # factor which differentiates this model by mclust
PXgood     <- array(0,c(n,k),dimnames=list(1:n,paste("group ",1:k,sep="")))
PXbad      <- array(0,c(n,k),dimnames=list(1:n,paste("group ",1:k,sep="")))

# ------------------------- #
# posteriors initialization #
# ------------------------- #

if(initialization=="random.soft"){
  
  z  <- array(runif(n*k),c(n,k)) # soft posterior probabilities (no-normalized) (n x k) 
  z  <- z/rowSums(z)             # soft posterior probabilities (n x k)
  for(j in 1:k){
    v[,2,j] <- runif(n,0,0.2)
    v[,1,j] <- 1-v[,2,j]
  }  
  
} 

if(initialization=="random.hard"){
  
  z  <- t(rmultinom(n, size = 1, prob=rep(1/k,k)))  # hard posterior probabilities (n x k)
  for(j in 1:k)
    v[,,j] <- t(rmultinom(n, size = 1, prob=c(0.9,0.1)))
  
} 

if(initialization=="manual"){ 
  
  z  <- start.z      # soft or hard classification
  v  <- start.v      # soft or hard classification of good and bad observations
  
}

if(initialization=="mclust"){
  
  mclustfit <- .gpcm(data=X, G=k, mnames=modelname, start=start, label=lab, pprogress=FALSE)
  for(j in 1:k){    
    mu[,j]        <- mclustfit$gpar[[j]]$mu
    Sigma[,,j]    <- mclustfit$gpar[[j]]$sigma
    invSigma[,,j] <- mclustfit$gpar[[j]]$invSigma
    #temp       <- eigen(Sigma[,,j])
    #Gamma[,,j] <- temp$vectors
    #Delta[,,j] <- diag(temp$values)/prod(temp$values)^(1/p)
    #lambda[j]  <- prod(temp$values)^(1/p)    
  }  
  z     <- mclustfit$z
  group <- apply(z,1,which.max)
  for(j in 1:k){ 
    Probj  <- (2*pi)^(-p/2)*det(Sigma[,,j])^(-1/2)*exp(-1/2*mahalanobis(x=X[which(group==j),], center=mu[,j], cov=invSigma[,,j], inverted=TRUE))
    #Probj  <- dmnorm(x=X[which(group==j),], mean = mu[,j], varcov=Sigma[,,j])
    wgood  <- Probj/max(Probj)
    v[which(group==j),1,j] <- (wgood)^0.2
    v[which(group==j),2,j] <- 1-(wgood)^0.2
  }
  
}

if(is.vector(label))
  z[ind.label,] <- class(label,k=k)

# ------------ #
# EM algorithm #
# ------------ #

# Preliminary definition of convergence criterions

check     <- 0
iteration <- 1
loglik    <- NULL
aloglik   <- NULL
aloglik   <- c(0,0)
a         <- NULL
a         <- c(0,0)

while(check<1){
  
  # ++++++ #
  # M-step #
  # ++++++ #
  
  # ------- #
  # Weights #
  # ------- #
  
  prior     <- colMeans(z)
  zv        <- z*v[,1,]
  zvcompl   <- z*v[,2,]
  if(!alphafix){
    if(alphacon){
      g <- function(alpha,j,z,v)      
        sum(z[,j]*(v[,1,j]*log(alpha)+v[,2,j]*log(1-alpha)))
      for(j in 1:k)
        priorgood[j] <- optimize(g,c(alphamin[j],1),maximum = TRUE,j=j,z=z,v=v)$maximum
      }
    if(!alphacon)
      priorgood <- colSums(zv)/colSums(z)
  }
  if(alphafix)
    priorgood <- alpha
  priorbad  <- 1-priorgood
    
  # ---------- #
  # mu & Sigma #
  # ---------- #

  correction <- v[,1,]+v[,2,]*matrix(rep(1/eta),n,k,byrow=TRUE)
  
  fitM <- .m.step(data=X, covtype=modelname, w=z, fact=correction, v=1, mtol=1e-10, mmax=10)
  for(j in 1:k){ 
    mu[,j]        <- fitM[[j]]$mu
    Sigma[,,j]    <- fitM[[j]]$sigma
    invSigma[,,j] <- fitM[[j]]$invSigma
    temp          <- eigen(Sigma[,,j])
    Gamma[,,j]    <- temp$vectors
    Delta[,,j]    <- diag(temp$values)/prod(temp$values)^(1/p)
    lambda[j]     <- prod(temp$values)^(1/p)
  }
  
  # ------------------- #
  # Inflation parameter #
  # ------------------- #  
  
  if(!etafix){
    if(etacon){
      f <- function(eta,j,p,zvcompl,X,mu,invSigma) 
        sum(-(p/2)*zvcompl[,j]*log(eta)-(1/2)*zvcompl[,j]*(1/eta)*mahalanobis(x=X, center=mu[,j], cov=invSigma[,,j], inverted=TRUE))
      for(j in 1:k)
        eta[j] <- optimize(f,c(1,etamax),maximum = TRUE,j=j,p=p,zvcompl=zvcompl,X=X,mu=mu,invSigma=invSigma)$maximum      
    }
    if(!etacon){
      for(j in 1:k)    
        eta[j] <- sum(zvcompl[,j]*mahalanobis(x=X, center=mu[,j], cov=invSigma[,,j], inverted=TRUE))/(p*sum(zvcompl[,j]))
    }
  }
    
  # ------- #
  # density #
  # ------- #
  
  zerostar <- 1e-323 # to avoid zero probabilities
  for(j in 1:k){
    PXgood[,j] <- (2*pi)^(-p/2)*(det(Sigma[,,j]))^(-1/2)*exp(-1/2*mahalanobis(x=X, center=mu[,j], cov=invSigma[,,j], inverted=TRUE)) 
    PXbad[,j]  <- (2*pi)^(-p/2)*(det(eta[j]*Sigma[,,j]))^(-1/2)*exp(-1/2*mahalanobis(x=X, center=mu[,j], cov=1/eta[j]*invSigma[,,j], inverted=TRUE))
    PXgood[,j] <- (PXgood[,j]<zerostar)*zerostar+(PXgood[,j]>=zerostar)*PXgood[,j]
    PXbad[,j]  <- (PXbad[,j]<zerostar)*zerostar+(PXbad[,j]>=zerostar)*PXbad[,j]
  }
   
  # ------------------------------------- # 
  # Global - Observed-data log-likelihood # 
  # ------------------------------------- #
  
  # model-based clustering
  
  if(!is.vector(label))
    llvalues <- sum(log( rowSums(matrix(rep(prior,n),n,k,byrow=TRUE)*(matrix(rep(priorgood,n),n,k,byrow=TRUE)*PXgood+matrix(rep(priorbad,n),n,k,byrow=TRUE)*PXbad))))
  
  # model-based classification
  
  if(is.vector(label)){    
    llvalueslab   <- z[ind.label,]*(log(matrix(rep(prior,nlab),nlab,k,byrow=T))+log( matrix(rep(priorgood,nlab),nlab,k,byrow=TRUE)*PXgood[ind.label,]+matrix(rep(priorbad,nlab),nlab,k,byrow=TRUE)*PXbad[ind.label,] ))
    llvaluesunlab <- log(rowSums(matrix(rep(prior,nunlab),nunlab,k,byrow=TRUE)*(matrix(rep(priorgood,nunlab),nunlab,k,byrow=TRUE)*PXgood[-ind.label,]+matrix(rep(priorbad,nunlab),nunlab,k,byrow=TRUE)*PXbad[-ind.label,])))
    llvalues      <- sum(llvalueslab)+sum(llvaluesunlab)    
  }
  
  loglik[iteration] <- llvalues
  
  # ----------------------------------------------- #
  # Aitkane's Acceleration-Based Stopping Criterion #
  # ----------------------------------------------- #
  
  if(iteration>2 & k > 1){
    if(abs(loglik[iteration-1]-loglik[iteration-2])>0){
      a[iteration-1]      <- (loglik[iteration]-loglik[iteration-1])/(loglik[iteration-1]-loglik[iteration-2])
      aloglik[iteration]  <- loglik[iteration-1]+(1/(1-a[iteration-1])*(loglik[iteration]-loglik[iteration-1]))
      if(abs(aloglik[iteration]-loglik[iteration])<threshold) 
        check <- 1
    }
    #if(abs(loglik[iteration-1]-loglik[iteration-2])==0)
    else
      check <- 1    
  }
  
  if(iteration==iter.max | k==1) check <- 1
  
  cat("*")
  #cat("Iteration",iteration,"\n")
  iteration <- iteration + 1
  
  # ++++++ #
  # E-Step #
  # ++++++ #
  
  z.num  <- matrix(rep(prior,n),n,k,byrow=TRUE)*(matrix(rep(priorgood,n),n,k,byrow=TRUE)*PXgood+matrix(rep(priorbad,n),n,k,byrow=TRUE)*PXbad)  # (n x k)
  z.den  <- rowSums(z.num)                        # n-vector
  z      <- z.num/matrix(rep(z.den,k),ncol=k)     # (n x k)
  
  if(is.vector(label))
    z[ind.label,] <- class(label,k=k)
    
  numgood <- matrix(rep(priorgood,n),n,k,byrow=TRUE)*PXgood
  numbad  <- matrix(rep(priorbad,n),n,k,byrow=TRUE)*PXbad
  v.den   <- matrix(rep(priorgood,n),n,k,byrow=TRUE)*PXgood+matrix(rep(priorbad,n),n,k,byrow=TRUE)*PXbad
  v[,1,]  <- numgood/v.den
  v[,2,]  <- numbad/v.den
  
}

  cat("\n")
  finalloglik <- loglik[iteration-1] 

# ---------------------------------------------------------------------------- #
# The EM-algorithm is finished                                                 #
# Check on the EM-monotonicity                                                 #
# Plot of the values of the observed-data log-likelihood versus the iterations #
# ---------------------------------------------------------------------------- #

# if(loglikplot){
#   
#   par(mai=c(0.84,0.8,0.012,0.004))
#   par(las = 3)
#   par(cex.axis=0.7)
#   par(cex.lab=1.2)
#   plot(0:(iteration-2),loglik[1:(iteration-1)],type="l",axes = FALSE,xlab="iterations",ylab="log-likelihood",lwd=2)
#   axis(1, at = 0:(iteration-2),label = 0:(iteration-2)) 
#   axis(2)
#   box(col = "black")
#   
# }

# --------------------- #
# Classification Matrix #
# --------------------- #

group <- apply(z,1,which.max)
innergroup  <- numeric(n)
for(i in 1:n)
  innergroup[i] <- ifelse(v[i,1,group[i]]<v[i,2,group[i]],"bad","*")
detection <- data.frame(group=group,innergroup=innergroup)

# -------------------- #
# Information criteria #
# -------------------- #
  
if(etafix){ 
  if(!alphafix)
    npar <- (k-1) + k + p*k + .ncovpar(modelname=modelname, p=p, G=k)
  if(alphafix)
    npar <- (k-1) + p*k + .ncovpar(modelname=modelname, p=p, G=k)  
}
if(!etafix){
  if(!alphafix)
    npar <- (k-1) + k + p*k + .ncovpar(modelname=modelname, p=p, G=k) + k
  if(alphafix)
    npar <- (k-1) + p*k + .ncovpar(modelname=modelname, p=p, G=k) + k
}

AIC  <- 2*finalloglik - npar*2
BIC  <- 2*finalloglik - npar*log(n)

z.const    <- (z<10^(-322))*10^(-322)+(z>10^(-322))*z   # vincolo per evitare i NaN nel calcolo di tau*log(tau)
hard.z     <- (matrix(rep(apply(z,1,max),k),n,k,byrow=F)==z)*1

if(is.vector(label)){  
  ECM       <- sum(hard.z[-ind.label,]*log(z.const[-ind.label,]))
  ICL       <- BIC+ECM  
}
if(!is.vector(label)){  
  ECM       <- sum(hard.z*log(z.const))
  ICL       <- BIC+ECM  
}
  
# ---- #
# plot #
# ---- #

# if(plot){
# 	if(loglikplot) 
#     x11()
# 	for(i in 2:p)
# 		for(j in 1:(i-1)){
# 			par(mai=c(0.84,0.8,0.012,0.004))
# 	  		par(las = 3)
# 	  		par(cex.axis=0.7)
# 	  		par(cex.lab=1.2)
# 	  		plot(X[,c(j,i)],col="white",xlab=paste("X_",j,sep=""),ylab=paste("X_",i,sep=""))
# 	  		text(X[,c(j,i)],labels=detection$innergroup,col=group,cex=0.7)
# 	  		box(col = "black")
# 			if((i+j)!=(2*p-1)) x11()
# 		}  
# }

result <- list(
  modelname = modelname,
  npar      = npar,
  X         = X,            
  k         = k,            
  p         = p,            
  n         = n,            
  prior     = prior,
  priorgood = priorgood,
  mu        = mu,
  Sigma     = Sigma,
  lambda    = lambda,
  Delta     = Delta,
  Gamma     = Gamma,
  eta       = eta,
  iter.stop = iteration,
  z         = z,
  v         = v,
  group     = group,
  detection = detection,
  loglik    = finalloglik,
  AIC       = AIC,
  BIC       = BIC,
  ICL       = ICL,          # alla McNicholas
  call      = match.call()
)

class(result) <- "MCM"
return(result)

alarm()

}

#####################
## model Selection ##
#####################

MS <- function(
  X,                            # matrix of data
  k,                            # vector of values for k
  model=NULL,                   # models to be considered in model selection
  initialization="mclust",      # initialization procedure: "random.soft", "random.hard", "manual", or "mclust"
  alphacon=TRUE,                # if TRUE, alpha is constrained to be >0.5
  alphamin=NULL,                # minimum proportion of good data in each group 
  alphafix=FALSE,               # if TRUE, the inner weights are estimated
  alpha=NULL,                   # vector of dimension k with proportion of good observations in each group
  etacon=TRUE,                  # if TRUE, eta is constrained to be >1
  etafix=FALSE,                 # if TRUE, eta is not estimated
  eta=NULL,                     # inflation parameters
  etamax=200,                   # maximum value of eta
  start.z=NULL,                 # (n x k)-matrix of soft or hard classification: it is used only if initialization="manual"  	
  start.v=NULL,                 # (n x 2 x k)-array of soft or hard classification in each group: it is used only if initialization="manual"  	
  start=0,                      # initialization for the package mixture
  ind.label=NULL,               # indexes of the labelled observations
  label=NULL,                   # groups of the labelled observations
  iter.max=1000,                # maximum number of iterations in the EM-algorithm
  threshold=1.0e-03             # stopping rule in the Aitken rule
)  
{
  call=match.call()

  if (is.data.frame(X)) 
    X <- as.matrix(X)  
  
  #   if(is.null(G)){
  #     G=1:3
  #   }
  #   if(is.null(modelNames)){
  #     modelnames=.cwm$modelNames   
  #   }
  #   if(is.null(method)){
  #     method="BIC"
  #   }
  
  n <- length(X)
  
  if(is.null(model))
    model <- c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVE","EEV","VVE","VEV","EVV","VVV")
  
  gridk     <- k
  numk      <- length(gridk)
  nummodel  <- length(model)
  
  IC <- array(0,c(numk,nummodel,4),dimnames=list(paste(gridk,"groups",sep=" "),paste("model",model,sep=" "),c("loglik","AIC","BIC","ICL")))
  
  cont <- 0
  par  <- list()
  #results <- NULL
  #pb  <- winProgressBar("Processing Simulations")
  #ee  <- proc.time()
  for(i in 1:numk){
    par[[i]] <- list()
    #results[[i]] <- NULL 
    for(j in 1:nummodel){
      cat("\n")
      cat(paste("Model ",model[j]," with ",gridk[i]," groups",sep=""))
      cat("\n")
      cont <- cont+1
      par[[i]][[j]] <- .MCM(
        X=X,			                      
        k=gridk[i],                            
        initialization=initialization,      
        modelname=model[j],                     
        alphacon=alphacon,                
        alphamin=rep(0.5,gridk[i]),
        alphafix=alphafix,               
        alpha=alpha,                   
        etacon=etacon,                  
        etafix=etafix,                 
        eta=eta,
        etamax=etamax,
        start.z=start.z,                 		
        start.v=start.v,                   	
        start=start,                      
        ind.label=ind.label,               
        label=label,                   
        iter.max=iter.max,                
        threshold=threshold           
        #loglikplot=FALSE,              
        #plot=FALSE                     
      )
      #results[[i]][[j]] <- temp
      IC[i,j,1] <- par[[i]][[j]]$loglik
      IC[i,j,2] <- par[[i]][[j]]$AIC
      IC[i,j,3] <- par[[i]][[j]]$BIC
      IC[i,j,4] <- par[[i]][[j]]$ICL
      
      #perc <- cont/(numk*nummodel)
      #setWinProgressBar(pb,perc)
      
    }
  }
  #close(pb)
  #proc.time()-ee
  
  cat("\n\n")
  cat("# ----------------------- #","\n")
  cat("# Model Selection Results #","\n")
  cat("# ----------------------- #","\n\n")
  
  # --- #
  # AIC #
  # --- #
  
  BestIndAIC <- which.max(IC[,,2])
  BestIndAIC <- arrayInd(BestIndAIC, .dim=c(numk,nummodel))
  BestAIC    <- IC[BestIndAIC[1],BestIndAIC[2],2]
  bestAIC    <- par[[BestIndAIC[1]]][[BestIndAIC[2]]]
  cat("Best AIC value of",BestAIC,"obtained for k =",gridk[BestIndAIC[1]],"group(s) with model",model[BestIndAIC[2]],"\n\n")
  
  # --- #
  # BIC #
  # --- #
  
  BestIndBIC <- which.max(IC[,,3])
  BestIndBIC <- arrayInd(BestIndBIC, .dim=c(numk,nummodel))
  BestBIC    <- IC[BestIndBIC[1],BestIndBIC[2],3]
  bestBIC    <- par[[BestIndBIC[1]]][[BestIndBIC[2]]]
  cat("Best BIC value of",BestBIC,"obtained for k =",gridk[BestIndBIC[1]],"group(s) with model",model[BestIndBIC[2]],"\n\n")
  
  # --- #
  # ICL #
  # --- #
  
  BestIndICL <- which.max(IC[,,4])
  BestIndICL <- arrayInd(BestIndICL, .dim=c(numk,nummodel))
  BestICL    <- IC[BestIndICL[1],BestIndICL[2],4]
  bestICL    <- par[[BestIndICL[1]]][[BestIndICL[2]]]
  cat("Best ICL value of",BestICL,"obtained for k =",gridk[BestIndICL[1]],"group(s) with model",model[BestIndICL[2]],"\n\n")
  
  bestk     <- bestmodel <- array(0,c(3),dimnames=list(c("AIC","BIC","ICL"))) 
  bestk     <- c(gridk[BestIndAIC[1]],gridk[BestIndBIC[1]],gridk[BestIndICL[1]]) 
  bestmodel <- c(model[BestIndAIC[2]],model[BestIndBIC[2]],model[BestIndICL[2]])
  best      <- data.frame(k=bestk,model=bestmodel)
  
  return(
    structure(
      list(
        call   = call,
        best   = best,
        #IC     = IC,
        bestAIC = bestAIC,
        bestBIC = bestBIC,
        bestICL = bestICL
        ),              
      class  = "MCMMS"
    )
  )  
  
}


