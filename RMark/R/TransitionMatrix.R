#' Multi-state Transition Functions
#' 
#' TransitionMatrix: Creates a transition matrix of movement parameters for a multi-state(strata)
#' model. It computes all Psi values for a multi-strata mark model and constructs a
#' transition matrix.  Standard errors and confidence intervals can also be
#' obtained.
#' 
#' find.possible.transitions:  Finds possible transitions; essentially it identifies where
#' stratum label A and B are in the same ch for all labels but the 
#' the transition could be from A to B or B to A or even ACB which is
#' really an A to C and C to B transition.
#' 
#' transition.pairs: Computes counts of transition pairs. The rows are the "from stratum" and the
#' columns are the "to stratum". So AB would be in the first row second column
#' and BA would be in the second row first column.  All intervening 0s are ignored.
#' These are transition pairs so AB0C is A to B and B to C but not A to C.
#' 
#' 
#' @param x Estimate table from \code{\link{get.real}} with a single record for
#' each possible transition
#' @param vcv.real optional variance-covariance matrix from the call to
#' \code{\link{get.real}}
#' @param ch vector of capture history strings for a multi-state analysis
#' @return TransitionMatrix: returns either a transition matrix (vcv.real=NULL) or a list of matrices
#' (vcv.real specified) named TransitionMat (transition matrix),
#' se.TransitionMat (se of each transition), lcl.TransitionMat (lower
#' confidence interval limit for each transition), and ucl.TransitionMat (upper
#' confidence interval limit for each transition). find.possible.transitions returns a 0/1 table where 1 means that t
#' both values are in one or more ch strings and transition.pairs returns a table of counts of transition pairs.
#' @author Jeff Laake
#' @usage  TransitionMatrix(x,vcv.real=NULL)
#' 
#'         find.possible.transitions(ch)
#' 
#'         transition.pairs(ch)
#' 
#' @aliases TransitionMatrix find.possible.transitions transition.pairs
#' @export TransitionMatrix find.possible.transitions transition.pairs
#' @seealso \code{\link{get.real}}
#' @keywords utility
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(mstrata)
#' # Show possible transitions in first 15 ch values
#' find.possible.transitions(mstrata$ch[1:15])
#' # Show transtion pairs for same data
#' transition.pairs(mstrata$ch[1:15])
#' #limit transtions to 2 and 3 character values for first 30 ch
#' transition.pairs(substr(mstrata$ch[1:30],2,3))
#' 
#' # fit the sequence of multistrata models as shown for ?mstrata
#' run.mstrata=function()
#' {
#' #
#' # Process data
#' #
#' mstrata.processed=process.data(mstrata,model="Multistrata")
#' #
#' # Create default design data
#' #
#' mstrata.ddl=make.design.data(mstrata.processed)
#' #
#' #  Define range of models for S; note that the betas will differ from the output
#' #  in MARK for the ~stratum = S(s) because the design matrix is defined using
#' #  treatment contrasts for factors so the intercept is stratum A and the other
#' #  two estimates represent the amount that survival for B abd C differ from A.
#' #  You can use force the approach used in MARK with the formula ~-1+stratum which
#' #  creates 3 separate Betas - one for A,B and C.
#' #
#' S.stratum=list(formula=~stratum)
#' S.stratumxtime=list(formula=~stratum*time)
#' #
#' #  Define range of models for p
#' #
#' p.stratum=list(formula=~stratum)
#' #
#' #  Define range of models for Psi; what is denoted as s for Psi
#' #  in the Mark example for Psi is accomplished by -1+stratum:tostratum which
#' #  nests tostratum within stratum.  Likewise, to get s*t as noted in MARK you
#' #  want ~-1+stratum:tostratum:time with time nested in tostratum nested in
#' #  stratum.
#' #
#' Psi.s=list(formula=~-1+stratum:tostratum)
#' #
#' # Create model list and run assortment of models
#' #
#' model.list=create.model.list("Multistrata")
#' #
#' # Add on a specific model that is paired with fixed p's to remove confounding
#' #
#' p.stratumxtime=list(formula=~stratum*time)
#' p.stratumxtime.fixed=list(formula=~stratum*time,fixed=list(time=4,value=1))
#' model.list=rbind(model.list,c(S="S.stratumxtime",p="p.stratumxtime.fixed",
#'                Psi="Psi.s"))
#' #
#' # Run the list of models
#' #
#' mstrata.results=mark.wrapper(model.list,data=mstrata.processed,ddl=mstrata.ddl)
#' #
#' # Return model table and list of models
#' #
#' return(mstrata.results)
#' }
#' mstrata.results=run.mstrata()
#' mstrata.results
#' # for the best model, get.real to get a list containing all Psi estimates
#' #  and the v-c matrix
#' Psilist=get.real(mstrata.results[[1]],"Psi",vcv=TRUE)
#' Psivalues=Psilist$estimates
#' # call Transition matrix using values from time==1; the call to the function
#' # must only contain one record for each possible transition. An error message is
#' # given if not the case
#' TransitionMatrix(Psivalues[Psivalues$time==1,])
#' # call it again but specify the vc matrix to get se and conf interval
#' TransitionMatrix(Psivalues[Psivalues$time==1,],vcv.real=Psilist$vcv.real)
#' }
TransitionMatrix=function(x,vcv.real=NULL)
{
#
#  Make sure certain conditions are met with argument values
#
  if(is.null(x$tostratum))
    stop("\nOnly works for Psi parameter of multistrata models\n")
  if(any(rowSums(table(list(x$stratum,x$tostratum)))!=(length(levels(x$stratum))-1)))
     stop("\nx should only contain one record for each stratum/tostratum pairing\n")
#
# Compute transition values for those estimated
#
  TransitionMat=tapply(x$estimate,list(x$stratum,x$tostratum),sum)
#
# next fill in the values computed by subtraction
#
  cc=col(TransitionMat)[is.na(TransitionMat)]
  rr=row(TransitionMat)[is.na(TransitionMat)]
  missing.index=cbind(rr,cc)[order(rr),]
  subtract.value=1-rowSums(TransitionMat,na.rm=TRUE)
  TransitionMat[missing.index]= subtract.value
#
# if the vcv.real matrix has been specified, then compute se and conf intervals
# for the transition matrix.
#
  if(!is.null(vcv.real))
  {
#
#    This code computes the se for the subtract stratum
#
     pin=tapply(x$par.index,list(x$stratum,x$tostratum),sum)
     se.TransitionMat=matrix(0,nrow=dim(pin)[1],ncol=dim(pin)[1])
     se.TransitionMat[!is.na(pin)]=x$se[match(pin[!is.na(pin)],x$par.index)]
     for (i in 1:dim(se.TransitionMat)[1])
     {
        vcv.indices=match(pin[i,][!is.na(pin[i,])],as.numeric(row.names(vcv.real)))
        se.TransitionMat[i,][is.na(pin[i,])]=sqrt(sum(vcv.real[vcv.indices,vcv.indices]))
     }
#
#    Next compute the conf interval for all the transitions using the
#    logit link individually for each transition.  This may not be
#    exactly correct but it matches what is done in MARK
#
     link.se=se.TransitionMat/(TransitionMat*(1-TransitionMat))
     link=log(TransitionMat/(1-TransitionMat))
     lcl.TransitionMat=inverse.link(link-1.96*link.se,"logit")
     ucl.TransitionMat=inverse.link(link+1.96*link.se,"logit")
#
#    Add strata.labels as row and column names for the matrices
#
     row.names(se.TransitionMat)=row.names(TransitionMat)
     colnames(se.TransitionMat)=row.names(TransitionMat)
     row.names(lcl.TransitionMat)=row.names(TransitionMat)
     colnames(lcl.TransitionMat)=row.names(TransitionMat)
     row.names(ucl.TransitionMat)=row.names(TransitionMat)
     colnames(ucl.TransitionMat)=row.names(TransitionMat)
     return(list(TransitionMat=TransitionMat,se.TransitionMat=se.TransitionMat,
         lcl.TransitionMat=lcl.TransitionMat,ucl.TransitionMat=ucl.TransitionMat))
  }
  else
     return(TransitionMat)
}

# Finds possible transitions; essentially it identifies where
# stratum label A and B are in the same ch for all labels but the 
# the transition could be from A to B or B to A or even ACB which is
# really an A to C and C to B transition.
find.possible.transitions=function(ch)
{
	splitch=sapply(ch,function(x) unlist(strsplit(x,split="")))
	strata=unique(as.vector(splitch))
	strata=strata[strata!="0"]
	strata=strata[order(strata)]
	trans=matrix(FALSE,nrow=length(strata),ncol=length(strata))
	i=0
	for(let1 in strata) 
	{ 
		i=i+1
		j=0
		for(let2 in strata)
		{
			j=j+1  
			trans[i,j]=as.numeric(length(grep(let2,ch[grep(let1,ch)]))>0)
		}
	}
	colnames(trans)=strata
	row.names(trans)=colnames(trans)
	trans[trans==0]=NA
	return(as.table(trans))
}
transition.pairs=function(ch)
{
	splitch=sapply(ch,function(x) unlist(strsplit(x,split="")))
	pairvals=function(x)
	{
		x=x[x!="0"]
		if(length(x)<=1)
			return(x)
		else
			return(paste(x[1:(length(x)-1)],x[2:length(x)],sep=""))
	}
	xx=unlist(apply(splitch,2,pairvals))
	if(!any(nchar(xx)>1))stop("\nNo transition pairs in data\n")
	xx=xx[nchar(xx)>1]
	xx=table(xx)
	let1=substr(names(xx),1,1)
	let2=substr(names(xx),2,2)
	xx=tapply(xx,list(let1,let2),mean)
	as.table(xx)
} 
