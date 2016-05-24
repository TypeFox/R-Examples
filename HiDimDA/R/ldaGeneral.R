### ldaGeneral.R  (2012-06-25)
###    
###
### Copyright 2012 A. Pedro Duarte Silva
###
### This file is part of the `HiDimDA' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

grpmeans <- function(x,grp) tapply(x,INDEX=grp,mean)
grpvar <- function(x,grp) tapply(x,INDEX=grp,var)
l2vnorm <- function(v)  sum(v^2)
normalize <- function(v,A)  v / sqrt(v%*%RightMult(A,v))

scalebygrps <- function(x,grouping,k=2,nk=NULL,n=NULL,p=NULL) 
{
	if (is.null(nk))  nk <- as.vector(table(grouping))
	if (is.null(n))  n <- sum(nk)
	if (is.null(p)) p <- ncol(x)
	vark <-  apply(x,2,grpvar,grp=grouping)
	globals <- sqrt(apply(matrix(rep(nk-1,p),k,p,byrow=FALSE)*vark,2,sum)/(n-k))
	list(Xscld=x/matrix(globals,n,p,byrow=TRUE),stdev=globals)
}

ldaGeneral <- function(data,grouping,prior,CorrAp,VSelfunct,ldafun,call,
			type=c("Dlda","Mlda","Slda","RFlda"),
#			PCAstep=NULL,q=NULL,maxcol=2000,nstarts=NULL,...)
			alwayscan,PCAstep=NULL,q=NULL,maxcol=2000,nstarts=NULL,...)
{
# Main routine to compute general linear discriminant analysis functions 

  type <- match.arg(type)

#  Compute the matrix Xw of within group deviations

  if  ( !is.matrix(data) || !is.factor(grouping)) stop("Arguments of wrong type.\n")
  if  ( type=="RFlda" && !is.numeric(q) ) 
	stop("When function ldaGeneral is called with type='RFlda' argument 'q' of  must be an integer.\n")
  grouping <- factor(grouping,levels=sort(unique(grouping))) 
  grpcds <- levels(grouping)

  nk <- table(grouping)
  n <- sum(nk)
  p <- ncol(data)
  k <- nrow(nk)
  df <- n-k
  if (min(nk)<=1) stop("There is (are) a group(s) with less than two observations.\n")
  if (n != length(grouping) || n != nrow(data) ) stop("Argument dimensions do not match.\n")
  if (is.null(colnames(data))) colnames(data) <- paste("X",1:p,sep="")
  if (prior[1] == "proportions") prior <- as.numeric(nk/n)
  names(prior) <- grpcds

  if (is.character(VSelfunct)) {
	if (VSelfunct=="none") {
		m <- p
		vkpt <- 1:p
	}
	else stop("Invalid value for VSelfunct argument.\n")
  }
  else  {
	if (!is.function(VSelfunct)) stop("Invalid value for VSelfunct argument.\n")
        SelV <- VSelfunct(data,grouping,...)
	if (!is.list(SelV)) stop("Invalid value for VSelfunct argument.\n")
	m <- SelV$nvkpt
	vkpt <- SelV$vkptInd
	data <- data[,vkpt,drop=FALSE]
  }  
  names(vkpt) <- colnames(data)

  if (CorrAp)  {
  	sclres <- scalebygrps(data,grouping,k=k,nk=nk,n=n,p=m)
  	data <- sclres$Xscld
  }
  if (type=="RFlda" && q > m) {
	warning("Number of assumed factors reduced from ",q," to the number of selected variables, that was equal to",m,"\n")
	q <- m
  }
  if (type=="Mlda")  {
	m0 <- m
 	if (PCAstep && m > n-1)  {
		uk0 <- apply(data,2,grpmeans,grp=grouping)
		FEgVcts <- rghtsngv(data/sqrt(n-1),nv=n-1)$v
		data <- data %*% FEgVcts   
		m <- n-1
	}
 	else if (!PCAstep && m > maxcol && m > n-1)  
		warning("Mlda will try to train the classifier with ",m," predictors. However,\n  with such a large number of variables it may run into memory problems.\n  Since there are ",n," observations and ",k," groups in the data, you can also\n  try to train the clasifier on its first ",n-k," Principal Components by seting\n  the 'PCAstep' argument to TRUE. That option is computationally cheaper\n  and usualy leads to comparable classification results.\n")
  }

  SigExpctRank <- min(m,df)
  uk <- apply(data,2,grpmeans,grp=grouping)
  dimnames(uk) <- list(levels(grouping),colnames(data))
  Xw <- matrix(0.,n,m)
  for (grp in 1:k)  {
	rowind <- which(grouping==grpcds[grp]) 
	Xw[rowind,] <- scale(data[rowind,],center=uk[grp,],scale=FALSE)
  }

#  Compute the matrices Sig and SigmaE of original and well-conditioned estimate total  within variances and covariances

   if (type=="Dlda") SigmaE <- DMat(colSums(Xw^2)/df)
   else  {
   	if ( !(type=="RFlda" && m>maxcol) && !(type=="Slda" && m>maxcol) ) Sig <- t(Xw) %*% Xw / (n - k)
   	if (type=="RFlda")  { 
		if (m>maxcol)  SigmaE <- FrobSigAp1(Xw/sqrt(df),min(n,m),q,nstarts)
   		else SigmaE <- FrobSigAp(Sig,q,nstarts)
	}	 
   }
   if (type=="Slda") {
	if (m>maxcol)  SigmaE <- ShrnkSigE(df,m,SigExpctRank,Sigma=NULL,SigmaSr=Xw/sqrt(df),check=FALSE,...)
	else SigmaE <- ShrnkSigE(df,m,SigExpctRank,Sigma=Sig,SigmaSr=NULL,check=FALSE,...)
   }
   if (type=="Mlda")  {
	if (ldafun=="canonical") {
		SigSigInvE <- MldaInvE(Sig,check=FALSE,onlyMinv=FALSE)
		SigmaE <- SigSigInvE$ME
		SigkptInv <- SigSigInvE$MInvE
	}
	else SigkptInv <- MldaInvE(Sig,check=FALSE)
   }
   else SigkptInv <- solve(SigmaE)

#  Compute the matrix Xdelta of between group deviations and the coefficients of the linear discriminant functions

   if (ldafun=="canonical")  {
  	HExpctRank <- min(m,k-1)
   	Xdelta <- matrix(0.,k,m)
   	for (grp in 1:k) Xdelta[grp,] <- uk[grp,] - colMeans(data)
#	if (m>maxcol)  {
#	}
#	else  {
		HSigkptInv <- RightMult(SigkptInv,t(Xdelta)%*%Xdelta/df)
		SepDec <- eigen(HSigkptInv)
		if (HExpctRank==1) Coef <- matrix(normalize(drop(Re(SepDec$vectors[,1])),A=SigmaE),m,1)
		else Coef <- apply(as.matrix(Re(SepDec$vectors[,1:HExpctRank])),2,normalize,A=SigmaE)
		sgvalues <- sqrt(Re(SepDec$values[1:HExpctRank]))
		names(sgvalues) <- colnames(Coef) <- paste("LD",1:HExpctRank,sep="")

#	}
   }
   else  {
 	Xavg <- matrix(0.,k-1,m)
   	Xdelta <- matrix(0.,k-1,m)
  	for (i in 2:k) {
		Xavg[i-1,] <- (uk[i,]+uk[1,])/2
   		Xdelta[i-1,] <- uk[i,] - uk[1,]
	}
   	Coef <- t(LeftMult(SigkptInv,Xdelta))
   	colnames(Coef) <- paste("CF",1:(k-1),sep="")
   	cnst <- array(dim=k-1)
   	for (grp in 2:k) cnst[grp-1] <-  matrix(Xavg[grp-1,],1,m) %*% Coef[,grp-1]
   } 
   if (type=="Mlda" && PCAstep && m0 > n-1)  {
	m <- m0
	uk <- uk0 
	Coef <- FEgVcts %*% Coef  
   }

   # return the results

   rownames(Coef) <- colnames(data)
   if (ldafun=="canonical")  {
   	if (type=="Dlda")  {
  		result <- list(prior=prior,means=uk,scaling=Coef,svd=sgvalues,vkpt=vkpt,nvkpt=m,N=n,call=call)
   		class(result) <- "canldaRes"
   	}   
   	else if (type=="Mlda")  {
   		if (CorrAp)  
     			result <- list(prior=prior,means=sweep(uk,2,sclres$stdev,FUN="*"),
					scaling=Coef/sclres$stdev,svd=sgvalues,vkpt=vkpt,nvkpt=m,N=n,call=call)
  		else result <- list(prior=prior,means=uk,scaling=Coef,svd=sgvalues,vkpt=vkpt,nvkpt=m,N=n,call=call)
   		class(result) <- "canldaRes"
   	}   
   	else if (type=="Slda")  {
   		if (CorrAp)  
     			result <- list(prior=prior,means=sweep(uk,2,sclres$stdev,FUN="*"),
					scaling=Coef/sclres$stdev,svd=sgvalues,
					vkpt=vkpt,nvkpt=m,SSig=SigmaE,SSigInv=SigkptInv,N=n,call=call)
  		else result <- list(prior=prior,means=uk,scaling=Coef,svd=sgvalues,
					vkpt=vkpt,nvkpt=m,SSig=SigmaE,SSigInv=SigkptInv,N=n,call=call)
   		class(result) <- c("Scanlda","canldaRes")
   	}   
   	else if (type=="RFlda")  {
   		if (CorrAp)  
     			result <- list(prior=prior,means=sweep(uk,2,sclres$stdev,FUN="*"),
					scaling=Coef/sclres$stdev,svd=sgvalues,vkpt=vkpt,nvkpt=m,
					q=q,SigFq=SigmaE,SigFqInv=SigkptInv,N=n,call=call)
  		else result <- list(prior=prior,means=uk,scaling=Coef,svd=sgvalues,vkpt=vkpt,nvkpt=m,
					q=q,SigFq=SigmaE,SigFqInv=SigkptInv,N=n,call=call)
   		class(result) <- c("RFcanlda","canldaRes")
  	}   
   }
   else if (ldafun=="classification")  {
   	if (type=="Dlda")  {
  		result <- list(prior=prior,means=uk,coef=Coef,cnst=cnst,vkpt=vkpt,nvkpt=m,N=n,call=call)
   		class(result) <- "clldaRes"
   	}   
   	else if (type=="Mlda")  {
   		if (CorrAp)  
     			result <- list(prior=prior,means=sweep(uk,2,sclres$stdev,FUN="*"),
					coef=Coef/sclres$stdev,cnst=cnst,vkpt=vkpt,nvkpt=m,N=n,call=call)
  		else result <- list(prior=prior,means=uk,coef=Coef,cnst=cnst,vkpt=vkpt,nvkpt=m,N=n,call=call)
   		class(result) <- "clldaRes"
   	}   
   	else if (type=="Slda")  {
   		if (CorrAp)  
     			result <- list(prior=prior,means=sweep(uk,2,sclres$stdev,FUN="*"),
					coef=Coef/sclres$stdev,cnst=cnst,
					vkpt=vkpt,nvkpt=m,SSig=SigmaE,SSigInv=SigkptInv,N=n,call=call)
  		else result <- list(prior=prior,means=uk,coef=Coef,cnst=cnst,
				vkpt=vkpt,nvkpt=m,SSig=SigmaE,SSigInv=SigkptInv,N=n,call=call)
   		class(result) <- c("Scllda","clldaRes")
   	}   
   	else if (type=="RFlda")  {
   		if (CorrAp)  
     			result <- list(prior=prior,means=sweep(uk,2,sclres$stdev,FUN="*"),
					coef=Coef/sclres$stdev,cnst=cnst,vkpt=vkpt,nvkpt=m,
					q=q,SigFq=SigmaE,SigFqInv=SigkptInv,N=n,call=call)
  		else result <- list(prior=prior,means=uk,coef=Coef,cnst=cnst,vkpt=vkpt,nvkpt=m,
					q=q,SigFq=SigmaE,SigFqInv=SigkptInv,N=n,call=call)
   		class(result) <- c("RFcllda","clldaRes")
	}
  }   
  result    # return(result)
}

LeftMult <- function(x,a) UseMethod("LeftMult") 

LeftMult.matrix <- function(x,a)
{
   p <- nrow(x)
   if ( (is.matrix(a) && p!=ncol(a)) || (!is.matrix(a) && p!=length(a)) ) stop("Argument dimensions do not match.\n")
   a %*% x  # return( a %*% x )
}

RightMult <- function(x,a) UseMethod("RightMult")  

RightMult.matrix <- function(x,a)
{
   p <- ncol(x)
   if ( (is.matrix(a) && p!=nrow(a)) || (!is.matrix(a) && p!=length(a)) ) stop("Argument dimensions do not match.\n")
   x %*% a  # return( x %*% a )
}

CovE <- function(object) UseMethod("CovE") 

CovE.Scanlda <- function(object) object$SSig
CovE.Scllda <- function(object) object$SSig
CovE.RFcanlda <- function(object) object$SigFq
CovE.RFcllda <- function(object) object$SigFq


ICovE <- function(object) UseMethod("ICovE") 

ICovE.Scanlda <- function(object) object$SSigInv
ICovE.Scllda <- function(object) object$SSigInv
ICovE.RFcanlda <- function(object) object$SigFqInv
ICovE.RFcllda <- function(object) object$SigFqInv

predict.canldaRes <- function(object,newdata,prior=object$prior,grpcodes=NULL,nbvrs=ncol(object$scaling),...)
{
   assigntomin <- function(scr) return(grpcodes[which.min(scr)])
   sqEuclidDist <- function(scr,grpmeans) apply(grpmeans,1,function(v1,v2) sum((v1-v2)^2),v2=scr) 	

   if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
   if (ncol(newdata)==1 && object$nvkp>1) newdata <- t(newdata)
   n <- nrow(newdata)
   k <- ncol(object$scaling) + 1
   if (is.null(grpcodes)) grpcodes <- 0:(k-1)
   if (is.null(rownames(newdata))) rownames(newdata) <- paste("Obs",1:n,sep="")
   if (is.null(rownames(newdata))) rownames(newdata) <- 1:n
   zmeans <- object$means%*%object$scaling[,1:nbvrs,drop=FALSE]
   if (is.null(object$vkpt)) zscores <- newdata%*%object$scaling[,1:nbvrs,drop=FALSE]
   else zscores <- newdata[,object$vkpt]%*%object$scaling[,1:nbvrs,drop=FALSE]
   sqrdist <- t(apply(zscores,1,sqEuclidDist,grpmeans=zmeans))
   dimnames(sqrdist) <- list(rownames(newdata),rownames(object$means))
   prioradj <- 2*log(prior)
   names(prioradj) <- rownames(object$means)
   res <- apply(scale(sqrdist,prioradj,scale=FALSE),1,assigntomin)   
   names(res) <- rownames(newdata)

   list(class=factor(res,levels=grpcodes),ZsqDistances=sqrdist,prior=prior,ZsqDprioradj=-prioradj,Z=zscores,Zmeans=zmeans)
}

predict.clldaRes <- function(object,newdata,prior=object$prior,grpcodes=NULL,...)
{
   assigntomax <- function(scr) return(grpcodes[which.max(scr)])

   if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
   if (ncol(newdata)==1 && object$nvkp>1) newdata <- t(newdata)
   n <- nrow(newdata)
   k <- ncol(object$coef) + 1
   if (is.null(grpcodes)) grpcodes <- 0:(k-1)
   if (is.null(rownames(newdata))) rownames(newdata) <- paste("Obs",1:n,sep="")
   cnst <- object$cnst + log(prior[1]/prior[-1])   # Adjust lda's treshold by the chosen (or estimated) prior probabilities
   if (k>2) {
 	if (is.null(object$vkpt)) scores <- cbind(rep(0.,n),scale(newdata%*%object$coef,center=cnst,scale=FALSE))
 	else scores <- cbind(rep(0.,n),scale(newdata[,object$vkpt]%*%object$coef,center=cnst,scale=FALSE))
	dimnames(scores) <- list(rownames(newdata),rownames(object$means))
  	res <- apply(scores,1,assigntomax)
   }
   else {
	if (is.null(object$vkpt)) scores <- scale(newdata%*%matrix(object$coef,object$nvkpt,1),center=cnst,scale=FALSE)
	else scores <- scale(newdata[,object$vkpt]%*%matrix(object$coef,object$nvkpt,1),center=cnst,scale=FALSE)
      	res <- rep(grpcodes[1],n)
	res[scores>0.] <- grpcodes[2]
	scores <- cbind(rep(0.,n),scores)
	names(scores) <- rownames(newdata)
	names(res) <- rownames(newdata)
   }
   list(class=factor(res,levels=grpcodes),x=scores[,-1])
}

print.canldaRes <- function(x,...)
{
	cat("Call:\n") ; print(x$call)  
	cat("\nPrior probabilities of groups:\n") ; print(x$prior)
	cat("\nGroup means:\n") ; print(x$means)
	cat("\nCoefficients of linear discriminants:\n") ; print(x$scaling) 
	if (length(x$prior)>2)  {
		egvalues <- x$svd^2	
		cat("Proportion of trace:\n")
		print(egvalues/sum(egvalues))
	}
	if (!(is.null(x$vkpt))) {
		cat("\nVariables kept in discriminant rule:\n") 
		if (!is.null(names(x$vkpt))) print(names(x$vkpt))
		else print(x$vkpt)
		cat("Number of variables kept:",x$nvkpt,"\n")
	}
}

print.clldaRes <- function(x,...)
{
	cat("Call:\n") ; print(x$call)  
	cat("\nPrior probabilities of groups:\n") ; print(x$prior)
	cat("\nGroup means:\n") ; print(x$means)
	cat("\nCoefficients of classification functions:\n") ; print(x$coef) 
	if (!(is.null(x$vkpt))) {
		cat("\nVariables kept in discriminant rule:\n") 
		if (!is.null(names(x$vkpt))) print(names(x$vkpt))
		else print(x$vkpt)
		cat("Number of variables kept:",x$nvkpt,"\n")
	}
}

coef.canldaRes <- function(object,...)  object$coef
coef.clldaRes <- function(object,...) object$coef



