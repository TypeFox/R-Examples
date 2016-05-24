### Dlda.R  (2012-06-23)
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

Dlda <- function(data,...) {
  if (is.null(class(data))) class(data) <- data.class(data)
  UseMethod("Dlda") 
}   

Dlda.default <- function(data,grouping,prior="proportions",VSelfunct=SelectV,ldafun=c("canonical","classification"),...)
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
  ldaGeneral(data,grouping,prior,FALSE,VSelfunct,type="Dlda",ldafun=ldafun,call=match.call(),...)  
}

Dlda.data.frame <- function(data,...)
{
   res <- Dlda.default(as.matrix(data),...)
   res$call <- match.call()
   res
}

is.Dlda <- function(x)  inherits(x,"Dlda")


