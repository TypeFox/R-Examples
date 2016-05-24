### Slda.R  (2012-06-19)
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

Slda <- function(data,...) {
  if (is.null(class(data))) class(data) <- data.class(data)
  UseMethod("Slda") 
}   

Slda.default <- function(data,grouping,prior="proportions",StddzData=TRUE,
			VSelfunct=SelectV,Trgt=c("CnstDiag","Idntty","VarDiag"),
			minp=20,ldafun=c("canonical","classification"),...)
{
  Trgt <- match.arg(Trgt)
  ldafun <- match.arg(ldafun)
  if (!is.matrix(data)) stop("'data' is not a matrix")
  if (!is.factor(grouping)) stop("'grouping' is not a factor")
  n <- nrow(data)
  if ( n != length(grouping)) stop("nrow(data) and length(grouping) are different")
  p <- ncol(data)
  k <- nrow(table(grouping))
  df <- n-k
  if (prior[1]!="proportions")  {
	if (!is.numeric(prior) || any(prior<0.||prior>1.) ) 
		stop("prior argument is not 'proportions' nor a vector of priors between 0. and 1.\n")
        lp <- length(prior) 
	if (lp!=k) stop(paste("Number of priors (",lp,") diferent than the number of groups (",k,").\n",sep=""))  
  }
  if (df<p) maxcol <- df
  else maxcol <- p	
  ldaGeneral(data,grouping,prior,StddzData,VSelfunct,type="Slda",ldafun=ldafun,call=match.call(),
		Trgt=Trgt,maxcol=maxcol,minp=minp,...)  
}

Slda.data.frame <- function(data,...)
{
   res <- Slda.default(as.matrix(data),...)
   res$call <- match.call()
   res
}

is.Slda <- function(x)  inherits(x,"Slda")

ShrnkSigE <- function(df,p,SigmaRank,Sigma=NULL,SigmaSr=NULL,check=TRUE,Trgt,
			minp=20,numtol=sqrt(.Machine$double.eps),...)
{
   if (p < minp && SigmaRank < p) stop("There are not enough variables to find reliable shrunk estimators of the covariance\n and not enough observations for non-singular unshrunken covariances.\n")
   if (is.null(Sigma) && is.null(SigmaSr) ) stop("Sigma and SigmaSr arguments of SrnkSigE1 cannot be simultaneously NULL.\n")
   if (!is.null(Sigma))  {
	if (check && !is.matrix(Sigma)) stop("SrnkSigE only accepts square matrix arguments.\n")
	if (check && p!=ncol(Sigma) || p!=nrow(Sigma)) stop("Wrong matrix dimensions.\n")
	if (check && SigmaRank > min(df,p)) 
		stop("Rank of covariance matrix cannot be above the number of available observations or variables.\n")
	if ( check && any(abs(Sigma-t(Sigma)) > numtol) ) stop("Original Sigma matrix is not symmetric.\n")
   	if (p<minp)  return(Sigma)	
	Sigmadec <- eigen(Sigma,symmetric=TRUE)
	D <- Sigmadec$values[1:SigmaRank]
	if (D[SigmaRank]<numtol) {
		D <- D[D>=numtol] 
		SigmaRank <- length(D)
	} 
	U <- Sigmadec$vectors[,1:SigmaRank]
   }
   else  {
   	if (p<minp)  return(SigmaSr %*% t(SigmaSr))
	SigSrsvddec <- rghtsngv(SigmaSr,nv=SigmaRank)
	D <- SigSrsvddec$d[1:SigmaRank]^2
	if (D[SigmaRank]<numtol) {
		D <- D[D>=numtol] 
		SigmaRank <- length(D)
	} 
	U <- SigSrsvddec$v[,1:SigmaRank]
   }
   trSigma <- sum(D)
   a1 <- trSigma/p	
   trSigma2 <- sum(D^2)
   a2 <- (df^2/((df-1)*(df+2)*p)) * (trSigma2-trSigma^2/df) 
   sqra1 <- a1^2 
   psqra1 <- p*sqra1 
   beta2 <- (a2+psqra1)/df
   if (Trgt[1]=="Idntty")  {
	delta2 <- ((df+1)*a2+psqra1)/df - 2*a1 + 1 
	lambda <- min(beta2/delta2,1.)  
   }
   if (Trgt[1]=="CnstDiag")  {
	delta2 <- ((df+1)*a2+psqra1)/df - sqra1
	lambda <- min(beta2/delta2,1.)  
	Trgt <- DMat(rep(a1,p))
   }
   if (Trgt[1]=="VarDiag")  {
	a2star <- (df/((df+2)*p)) * sum(diag(Sigma)^2)
	eta2 <- -2*a2star/df
	delta2 <- ((df+1)*a2+psqra1-(df+2)*a2star)/df 
	lambda <- min((beta2+eta2)/delta2,1.)  
	if (!is.null(Sigma)) Trgt <- DMat(diag(Sigma))
	else Trgt <- DMat(apply(SigmaSr,2,l2vnorm))
    }
    if (lambda==1.)  {
	if (Trgt[1]=="Idntty") return(DMat(rep(1.,p)))
	else return(Trgt)
    }
    ShrnkMat(Trgt=Trgt,U=U,D=D,p=p,q=SigmaRank,Intst=lambda)
}

