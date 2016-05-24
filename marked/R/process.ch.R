#' Process release-recapture history data
#' 
#' Creates needed constructs from the release-recapture history.
#' 
#' 
#' @param ch Vector of character strings; each character string is 
#' composed of either a constant length sequence of single characters (01001) or 
#' the character string can be comma separated if more than a single 
#' character is used (1S,1W,0,2W). If comma separated, each string must contain a constant
#' number of elements. 
#' @param freq Optional vector of frequencies for ch; if missing assumed to be
#' a; if <0 indicates a loss on capture
#' @param all FALSE is okay for cjs unless R code used to compute lnl instead of
#' FORTRAN; must be true for js because it returns additional quantities needed for entry prob.
#' @return \item{nocc}{number of capture occasions} \item{freq}{absolute value
#' of frequency for each ch} \item{first}{vector of occasion numbers for first
#' 1} \item{last}{vector of occasion numbers for last 1} \item{loc}{vector of
#' indicators of a loss on capture if set to 1} \item{chmat}{capture history
#' matrix} \item{FtoL}{1's from first (1) to last (1) and 0's elsewhere; only
#' if all==TRUE} \item{Fplus}{1's from occasion after first (1) to nocc(last
#' occasion); only if all==TRUE} \item{Lplus}{1's from occasion after last (1)
#' to nocc; only if all==TRUE} \item{L}{1's from last (1) to nocc; only if
#' all==TRUE}
#' @export
#' @author Jeff Laake
process.ch=function(ch,freq=NULL,all=FALSE)
#################################################################################
# process.ch - from a capture history creates vector of first and last times seen,
#              freq vector and indicator matrices used in log-likelihood calculations
#
#  Argument:
#        ch       - vector of character strings of 0/1
#        freq     - frequency of that ch; if null assumed to be 1; if <0 it
#                   signifies a loss on capture at last capture event
#        all      - if TRUE, computes all indicator matrices
#
#  Value: list with following elements
#        nocc        - number of capture occasions
#        freq        - absolute value of frequency for each ch
#        first       - vector of occasion numbers for first 1
#        last        - vector of occasion numbers for last 1
#        loc         - indicator of a loss on capture if set to 1
#        chmat       - capture history matrix
#        The following only returned if all==TRUE
#        FtoL        - 1's from first (1) to last (1) and 0's elsewhere (excluding
#        Fplus       - 1's from occasion after first (1) to nocc(last occasion)
#        Lplus       - 1's from occasion after last (1) to nocc
#        L           - 1's from last (1) to nocc
#        First       - 1's from occasion first (1) to nocc(last occasion)
#################################################################################
{
#  is ch comma separated? If not, separate by commas
   if(length(grep(",",ch[1]))==0)
		ch=sapply(strsplit(ch,""),paste,collapse=",")
   ch.lengths=sapply(strsplit(ch,","),length)
   nocc=ch.lengths[1]
   if(any(ch.lengths!=nocc))
	   stop("\nCapture history length is not constant. \nch must be a character string with constant length or comma separated with constant number of elements \n")	
   nch=length(ch)
   if(is.null(freq))freq=rep(1,nch)
# in case multistate data are passed change all non-zero to 1
   chmat=matrix((unlist(strsplit(ch,","))),byrow=TRUE,ncol=nocc,nrow=nch)
   ch=apply(t(apply(splitCH(ch),1,function(x){ 
							   x[x!="0"]=1 
							   return(x)
						   })),1,paste,collapse="")
#  create a matrix with 1:nocc in each row and one row for each ch
   nums=matrix(1:nocc,nrow=nch,ncol=nocc,byrow=TRUE)
#  store in a temp matrix and assign any 0 value to NA
   ymat=matrix(as.numeric(unlist(strsplit(ch,""))),byrow=TRUE,ncol=nocc,nrow=nch)
   if(suppressWarnings(all(is.numeric(as.numeric(chmat)))))chmat=ymat
   ymat[ymat==0]=NA
#  multiply nums matrix times the chmat
   nums=nums*ymat
#  use apply to get the minimum occasion and max occasion excluding NA values
   first=apply(nums,1,min,na.rm=TRUE)
   last=apply(nums,1,max,na.rm=TRUE)
   loc=rep(0,length(first))
   loc[freq<0]=1
   freq=abs(freq)
#  using the first and last values for each ch, compute the indicator matrices as defined above
   if(all)
   {
      FtoL=t(apply(cbind(first,last),1,function(x,nocc) {return(c(rep(0,x[1]),rep(1,x[2]-x[1]),rep(0,nocc-x[2])))},nocc=nocc))
      First=t(sapply(first,function(x,nocc){return(c(rep(0,x[1]-1),rep(1,nocc-x[1]+1)))},nocc=nocc))
      Fplus=t(sapply(first,function(x,nocc){return(c(rep(0,x[1]),rep(1,nocc-x[1])))},nocc=nocc))
      Lplus=t(sapply(last,function(x,nocc){return(c(rep(0,x[1]),rep(1,nocc-x[1])))},nocc=nocc))
      L=t(sapply(last,function(x,nocc){return(c(rep(0,x[1]-1),rep(1,nocc-x[1]+1)))},nocc=nocc))
#  return a list with each of the values
      return(list(nocc=nocc,freq=freq,first=first,last=last,loc=loc,chmat=chmat,FtoL=FtoL,Fplus=Fplus,Lplus=Lplus,L=L,First=First))
   }
   else
      return(list(nocc=nocc,freq=freq,first=first,last=last,loc=loc,chmat=chmat))      
}
