### Mlda.R  (2012-06-23)
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

Mlda <- function(data,...) {
  if (is.null(class(data))) class(data) <- data.class(data)
  UseMethod("Mlda") 
}   

Mlda.default <- function(data,grouping,prior="proportions",StddzData=TRUE,VSelfunct=SelectV,
ldafun=c("canonical","classification"),PCAstep=FALSE,...)
{
  ldafun <- match.arg(ldafun)
  if (!is.matrix(data)) stop("'data' is not a matrix")
  if (!is.factor(grouping)) stop("'grouping' is not a factor")
  n <- nrow(data)
  if ( n != length(grouping)) stop("nrow(data) and length(grouping) are different")
  if (prior[1]!="proportions")  {
	if (!is.numeric(prior) || any(prior<0.||prior>1.) ) 
		stop("prior argument is not 'proportions' nor a vector of priors between 0 and 1.\n")
	k <- nrow(table(grouping))
        lp <- length(prior) 
	if (lp!=k) stop(paste("Number of priors (",lp,") diferent than the number of groups (",k,").\n",sep=""))  
  }
  ldaGeneral(data,grouping,prior,StddzData,VSelfunct,type="Mlda",ldafun=ldafun,call=match.call(),PCAstep=PCAstep,...)  
}

Mlda.data.frame <- function(data,...)
{
   res <- Mlda.default(as.matrix(data),...)
   res$call <- match.call()
   res
}

is.Mlda <- function(x)  inherits(x,"Mlda")

MldaInvE <- function(M,check=TRUE,onlyMinv=TRUE,numtol=sqrt(.Machine$double.eps))
{
	if (check && !is.matrix(M)) stop("NldaInvE only accepts square matrix arguments.\n")
	p <- ncol(M)
	if (check && p!=nrow(M)) stop("NldaInvE only accepts square matrix arguments.\n")
	if ( check && any(abs(M-t(M)) > numtol) ) stop("Original M matrix is not symmetric.\n")
   	SpeDec <- eigen(M,symmetric=TRUE)
	if (check && SpeDec$values[p] < -numtol) 
		stop("NldaInvE only accepts positive definite or positive semidefinite matrix arguments.\n")
        muegval <- mean(SpeDec$values)
        newegval <- sapply(SpeDec$values,function(x,y) return(max(x,y)),y=muegval)
	MInvE <- matrix(0.,p,p)
	if (!onlyMinv)  ME <- matrix(0.,p,p)
	for (i in 1:p) {
		vivit <- outer(SpeDec$vectors[,i],SpeDec$vectors[,i]) 
		MInvE <- MInvE + vivit/newegval[i] 
		if (!onlyMinv) ME <- ME + newegval[i]*vivit
	} 
	if (!onlyMinv)  return(list(ME=ME,MInvE=MInvE))
	MInvE   # return(MInvE)
}

