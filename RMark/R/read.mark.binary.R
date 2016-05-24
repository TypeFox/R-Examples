#' Reads binary file output from MARK and returns a list of the results
#' 
#' Window and linux versions to read binary files created by MARK
#' 
#' 
#' @aliases read.mark.binary read.mark.binary.linux
#' @param filespec Filename specification for binary output file from
#' MARK;named here as markxxx.vcv
#' @param derived_labels vector of labels for derived parameters; NULL if no derived parameters for model
#' @return List of estimates, se, lcl, ucl and var-cov matrices for beta, real
#' and derived estimates \item{beta}{Dataframe for beta parameters containing
#' estimates, se, lcl, ucl} \item{beta.vcv}{variance-covariance matrix for beta
#' estimates} \item{real}{Dataframe for real parameters containing estimates,
#' se, lcl, ucl} \item{real.vcv}{variance-covariance matrix for real estimates}
#' \item{derived}{Dataframe for derived parameters (if any) containing
#' estimates, se, lcl, ucl} \item{derived.vcv}{variance-covariance matrix for
#' derived estimates (if any)}
#' @author Jeff Laake
#' @seealso \code{\link{extract.mark.output}}
#' @keywords utility
read.mark.binary <-
function(filespec,derived_labels)
{
z=file(filespec,"rb")
beta=list()
whatread=readChar(z,12)
ncovs=readBin(z,integer(),1)
nlogit=readBin(z,integer(),1)
ngrps=readBin(z,integer(),1)
nests=readBin(z,integer(),1)
#whatread=readBin(z,integer(),floor(max(ncovs, nlogit, ngrps, nests)/2)*2)
#whatread=readChar(z,8)
whatread=""
while(whatread!="B")
  whatread=readChar(z,1)
#if(whatread!="BETA PAR")
whatread=readChar(z,7)
beta$estimate=readBin(z,numeric(),ncovs)
whatread=readChar(z,16)
beta$se=readBin(z,numeric(),ncovs)
whatread=readChar(z,16)
beta$lcl=readBin(z,numeric(),ncovs)
whatread=readChar(z,16)
beta$ucl=readBin(z,numeric(),ncovs)
beta.vcv=matrix(0,nrow=ncovs,ncol=ncovs)
for (i in 1:ncovs)
{
   whatread=readChar(z,16)
   beta.vcv[i,]=readBin(z,numeric(),ncovs)
}
whatread=readChar(z,16)
real=list()
real$estimate=readBin(z,numeric(),nlogit)
whatread=readChar(z,16)
real$se=readBin(z,numeric(),nlogit)
whatread=readChar(z,16)
real$lcl=readBin(z,numeric(),nlogit)
whatread=readChar(z,16)
real$ucl=readBin(z,numeric(),nlogit)
real.vcv=matrix(0,nrow=nlogit,ncol=nlogit)
for (i in 1:nlogit)
{
whatread=readChar(z,16)
real.vcv[i,]=readBin(z,numeric(),nlogit)
}
if(!is.null(derived_labels))
{
	derived=vector("list",length(derived_labels))
	names(derived)=derived_labels
	derived.vcv=vector("list",length(derived_labels))
	names(derived.vcv)=derived_labels
	for(label in derived_labels)
	{
		whatread=readChar(z,16)
		nderiv=readBin(z,integer(),1)
		if(length(nderiv)!=0)
		{
			derived[[label]]=list()  
			whatread=readChar(z,16)
			derived[[label]]$estimate=readBin(z,numeric(),nderiv)
			whatread=readChar(z,16)
			derived[[label]]$se=readBin(z,numeric(),nderiv)
			whatread=readChar(z,16)
			derived[[label]]$lcl=readBin(z,numeric(),nderiv)
			whatread=readChar(z,16)
			derived[[label]]$ucl=readBin(z,numeric(),nderiv)
			derived[[label]]=as.data.frame(derived[[label]])
			derived.vcv[[label]]=matrix(0,nrow=nderiv,ncol=nderiv)
			for (i in 1:nderiv)
			{
				whatread=readChar(z,16)
				value=readBin(z,numeric(),nderiv)
				if(length(value)!=0) derived.vcv[[label]][i,]=value
			}
		}
	}
} else
{
 derived=NULL
 derived.vcv=NULL
}
close(z)
return(list(beta=as.data.frame(beta),beta.vcv=beta.vcv,real=as.data.frame(real),real.vcv=real.vcv,derived=derived,derived.vcv=derived.vcv))
}
