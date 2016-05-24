#' Multivariate Multistate (mvms) Specification
#' 
#' Creates list data structure from mvms specification 
#' 
#' Accepts a mvms specification which is a list with named character vectors and
#' optionally a vector named exclude. The length of the list (except for exclude)
#' is the multivariate dimension. The name of each dimension is the list element name.  
#' Each character vector specifies the one character labels for the states of that dimension and optionally a 
#' reserved character "u" to specify that there is state uncertainty for that
#' dimension. The name vector exclude can be used to remove variable combinations that cannot 
#' occur in the data.
#' 
#' The code tests to make sure that the input mvms specification is of the correct 
#' structure and it stops with an error message if not. The code returns a list
#' structure with a number of elements described under return value below.
#' 
#' @param x a multivariate multistate (mvms) specification as described above
#' @return a list with the following elements: 1) mvms - the input specification,
#' 2) nd - the number of dimensions, 3) df - the dataframe containing all combinations 
#' of observations across dimensions including uncertain states, 4) df.states - the dataframe with all
#' combinations of states across dimensions, 5) uncertain - boolean vector with nd elements
#' indicating whether there is uncertainy in states for each dimension.
#' @author Jeff Laake
#' @export 
#' @examples
#' set_mvms(list(location=c("A","B","C"),repro_status=c("N","P","u"),exclude=c("CP")))
set_mvms=function(x)
{
	if(!is.list(x))
		stop("mvms specification must be a list")
	if(any(names(x)==""))
	    stop("all dimensions of mvms specification must be named")
	if(any(!sapply(x,function(x) is.vector(x))))
		stop("each set of states must be specified as a vector")
	if(any(!sapply(x,function(x) mode(x)=="character")))
	    stop("each set of states must be specified as a character vector")
	if(!is.null(x$exclude))
	{
		exclude=x$exclude
		x$exclude=NULL
	}
	else
		exclude=NULL
	x=lapply(x,function(x)x[order(x)])
	nd=length(x)
	uncertain=sapply(x,function(x)"u"%in%x)
	nou=lapply(x,function(x) x[x!="u"])
	df=expand.grid(x)
	strata=apply(df,1,paste,collapse="")
	df=df[order(strata),,drop=FALSE]
	df.nou=expand.grid(nou)
	if(!is.null(exclude))
	{
		df=df[!strata%in%exclude,,drop=FALSE]
		strata=apply(df.nou,1,paste,collapse="")
		df.nou=df.nou[!strata%in%exclude,,drop=FALSE]
	}
	strata=apply(df.nou,1,paste,collapse="")
	df.nou=df.nou[order(strata),,drop=FALSE]
	return(list(mvms=x,nd=nd,df=df,df.states=df.nou,uncertain=uncertain))
}
#' Multivariate Multistate (mvms) Design Data 
#' 
#' Creates a dataframe with design data for MvMS model for a single occasion 
#' @param df.states dataframe of states created from set_mvms
#' @param df is dataframe of observations; provided for parameter delta in which both state and observation are needed in design data
#' @param transition if TRUE, creates design data for a state transition (from to); otherwise just state variables
#' @return a dataframe to be used with design data for mvms model
#' @author Jeff Laake
#' @export 
#' @examples
#' x=set_mvms(list(location=c("A","B","C"),repro_status=c("N","P","u")))
#' mvms_design_data(x$df.states)
#' mvms_design_data(x$df.states,transition=FALSE)
#' 
mvms_design_data=function(df.states,df=NULL,transition=TRUE)
{
	z=apply(df.states,1,paste,collapse="")
	if(transition)
	{
		z=expand.grid(z,z)
		colnames(z)=c("stratum","tostratum")
	    for(i in 1:2)
	    {
		    x=as.data.frame(do.call("rbind",strsplit(t(sapply(z[,i],as.character)),"")))
			if(ncol(x)>1)
			{
				if(i==1)
					colnames(x)=names(df.states)
				else
					colnames(x)=paste("to",names(df.states),sep="")
				z=cbind(z,x)
			}
	    }
	    return(z)
    } else
	{
		if(is.null(df))
		{
			if(ncol(df.states)>1)
				return(cbind(data.frame(stratum=z),df.states))
			else
				return(cbind(data.frame(stratum=z)))
		} else
		{
			uncertain=sapply(df,function(x) length(grep("u",x))>0)
			if(ncol(df.states)>1)
				xx=cbind(data.frame(stratum=z),df.states)
			else
				xx=cbind(data.frame(stratum=z))
			assignobs=function(x)
			{
				xx=as.list(as.data.frame(rbind(as.matrix(x)),stringsAsFactors=FALSE))
				xx=c(x[!uncertain],lapply(x[uncertain],function(x) c(x,"u")))
				xx=expand.grid(xx)
				xx=xx[order(apply(xx,1,paste,collapse="")),,drop=FALSE]
			}
			if(ncol(df.states)>1)
			    obs=do.call("rbind",apply(xx[,-1],1, assignobs))
			else
				obs=do.call("rbind",apply(xx,1, assignobs))
	        nobs=nrow(obs)/nrow(xx)
			colnames(obs)=paste("obs",names(obs),sep=".")
			return(cbind(xx[rep(1:nrow(xx),each=nobs),,drop=FALSE],obs))
		}
	}
}

