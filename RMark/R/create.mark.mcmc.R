#' Create mcmc object for analysis with coda
#' 
#' Reads in mcmc file from program MARK in binary and returns an mcmc object that can be used with coda functions which are 
#' most easily accessed via codamenu().
#' 
#' @param filename name of file containing mcmc values
#' @param ncovs number of covariates
#' @param nmeans number of means
#' @param ndesigns number of designs
#' @param nsigmas number of sigmas
#' @param nrhos number of rhos
#' @param nlogit number of logits
#' @param include if TRUE it includes ir/propJumps fields
#' @return An mcmc object if one chain and an mcmc list if more than one chain. Each can
#' be used with the coda package
#' @author Jeff Laake
#' @export
#' @import coda
#' @keywords utility
create.mark.mcmc <-function(filename,ncovs,nmeans,ndesigns,nsigmas,nrhos,nlogit,include=F)
{
#
#  Compute nestimates and the number of values per record (nvals)
#
nestimates=(ncovs+nmeans+ndesigns+nsigmas+nrhos);
nvals=4+2*nestimates+nlogit;

#
# create column names depending on values of arguments
#
   cnames=c("PropJumps",paste("ir",1:nestimates,sep=""),paste("Beta",1:ncovs,sep=""))
   if(nmeans>0)cnames=c(cnames,paste("mean",1:nmeans,sep=""))
   if(ndesigns>0)cnames=c(cnames,paste("design",1:ndesigns,sep=""))
   if(nsigmas>0)cnames=c(cnames,paste("sigma",1:nsigmas,sep=""))
   if(nrhos>0)cnames=c(cnames,paste("rho",1:nrhos,sep=""))
   cnames=c(cnames,paste("M2logL",sep=""))
   cnames=c(cnames,paste("real",1:nlogit,sep=""))

#
# Open file
#
z=file(filename,"rb")
xmcmc=matrix(NA,nrow=1000,ncol=nvals)
size=1000
#
# While the number of records read is less than 1000, keep reading in data and appending to mcmc matrix
#
# 
i=1
while(1==1)
{
   for(j in 1:4) dummy=readBin(z,raw())
   if(length(dummy)==0)break
   xmcmc[i,]=readBin(z,numeric(),n=nvals)
   for(j in 1:4) dummy=readBin(z,raw())
   i=i+1 
   if(i>size)
   {
      xmcmc=rbind(xmcmc,matrix(NA,nrow=1000,ncol=nvals))
	  size=size+1000
   }
}
xmcmc=xmcmc[1:(i-1),]
#
# See if there is more than one MC chain; if so create a *list* of mcmc objects
#

if(max(xmcmc[,1])!=min(xmcmc[,1]))
{
   mcmclist=list()
   for(i in 1:max(xmcmc[,1]))
   {
       if(include) {
          xmat=xmcmc[xmcmc[,1]==i,-(1:2)] 
          }
       else
       {
          xmat=xmcmc[xmcmc[,1]==i,-(1:(nestimates+3))]
       }
       # now assign column names
       if(include) {
       colnames(xmat)=cnames 
       }
       else
       { colnames(xmat)=cnames[-(1:(nestimates+1))] }
       mcmclist[[i]]=coda::mcmc(xmat)
   }
   mcmclist=coda::mcmc.list(mcmclist)
}
else
#
# MCMC input is only a single chain; add names here and return value
#
{
   if(include) {
      xmat=xmcmc[,-(1:2)]
      }
   else
   {
      
      xmat=xmcmc[,-(1:(nestimates+3))]
      cnames=cnames[-(1:(nestimates+1))]
   }

   colnames(xmat)=cnames
   mcmclist=coda::mcmc(xmat)
}
#
# return result as an mcmc object or list dataframe
#
message("\n");
message("+----------------------------------------------------------------+\n");
message("| Conversion of MARK MCMC output for CODA processing complete... |\n");
message("+----------------------------------------------------------------+\n\n");

message("(The mcmc.bin.read.R script has created a MCMC object called \n");
message(" 'mcmcdata' which can now be processed by the CODA package...if new\n");
message("  to R, this is most easily accomplished by invoking the codamenu() \n")
message("  command from within the R console.\n\n");

return(mcmclist)

} 
                                               
