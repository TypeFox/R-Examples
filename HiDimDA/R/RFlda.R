### RFlda.R  (2012-06-23)
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

RFlda <- function(data,...) {
  if (is.null(class(data))) class(data) <- data.class(data)
  UseMethod("RFlda") 
}   

RFlda.default <- function(data,grouping,q=1,prior="proportions",CorrAp=TRUE,maxq=5,VSelfunct=SelectV,
			ldafun=c("canonical","classification"),
			nstarts=1,CVqtrials=1:3,CVqfolds=3,CVqrep=1,CVqStrt=TRUE,...)
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
  maxcol <- 2000   # Maximum dimensionality allowed for square matrices. Matrix square roots are ued when this limit is surpassed.
  if (is.numeric(q)) {
	if (q > maxq) stop("The number of factors exceeds its upper limit of",maxq," . This limit can be increased by changing the value of the 'maxq' argument, but the resulting classification rule may take too long to compute.")
	if (q < 0) stop("Argument 'q' must be a non-negative integer\n.")
	if (q==0) {
		warning("Since argument 'q' was set to 0 a diagonal linear discriminant rule will be assumed\n")
		return(ldaGeneral(data,grouping,prior,FALSE,VSelfunct,type="Dlda",ldafun=ldafun,call=match.call()))
	}
	return(ldaGeneral(data,grouping,prior,CorrAp,VSelfunct,type="RFlda",ldafun=ldafun,call=match.call(),
			q=q,maxcol=maxcol,nstarts=nstarts,...))
  }
  else  {

  	TAl <- function(clrule,...) clrule

	if (q!="CVq") stop("Argument q must be numeric or equal to the string 'CVq'")
	if (max(CVqtrials) > maxq) stop("The maximum number of factors tp be tested exceeds the upper limit of",maxq," . This limit can be increased by changing the value of the 'maxq' argument, but the resulting classification rule may take too long to compute.")
  	nbqtrials <- length(CVqtrials)
  	Results <- vector("list",nbqtrials)
  	CVResults <- array(dim=nbqtrials)
  	totrep <- CVqrep*CVqfolds
  	gcds <- levels(grouping)
  	for (q in CVqtrials)  Results[[q]] <- ldaGeneral(data,grouping,prior,CorrAp,VSelfunct,type="RFlda",ldafun=ldafun,call=match.call(),
								q=q,maxcol=maxcol,nstarts=nstarts,...)
  	for (q in CVqtrials)  {   
		tmp <- DACrossVal(data,grouping,TrainAlg=TAl,kfold=CVqfolds,CVrep=CVqrep,Strfolds=CVqStrt,
					ldafun="classification",call=match.call(),grpcodes=gcds,clrule=Results[[q]],...)
		errates <- apply(tmp[,,1]*tmp[,,2],2,sum)/n
		CVResults[q] <- mean(errates)
    	}
    	bestq <- which.min(CVResults)	
    	Results[[bestq]]  # return(Results[[bestq]])
    }
}

RFlda.data.frame <- function(data,...)
{
   res <- RFlda.default(as.matrix(data),...)
   res$call <- match.call()
   res
}

is.RFlda <- function(x)  inherits(x,"RFlda")


