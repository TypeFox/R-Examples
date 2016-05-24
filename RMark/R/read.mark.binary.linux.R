read.mark.binary.linux <-
function(filespec,derived_labels)
#
# read.mark.binary 
#
# Arguments:
# 
# Value: list of estimates, se, lcl, ucl and var-cov matrices
#        for beta, real and derived estimates
{
z=file(filespec,"rb")
beta=list()
whatread=readBin(z,raw(),12)
ncovs=readBin(z,integer(),1,size=4)
nlogit=readBin(z,integer(),1,size=4)
ngrps=readBin(z,integer(),1,size=4)
nests=readBin(z,integer(),1,size=4)
whatread=""
while(whatread!="B" & whatread!="E")
{
  whatread=try(readChar(z,1),silent=TRUE)
  if(class(whatread)=="try-error") whatread=""
}
if(whatread=="B")
  whatread=readBin(z,raw(),7)
else
	whatread=readBin(z,raw(),6)
beta$estimate=readBin(z,numeric(),ncovs,size=8)
whatread=readBin(z,raw(),16)
beta$se=readBin(z,numeric(),ncovs,size=8)
whatread=readBin(z,raw(),16)
beta$lcl=readBin(z,numeric(),ncovs,size=8)
whatread=readBin(z,raw(),16)
beta$ucl=readBin(z,numeric(),ncovs,size=8)
beta.vcv=matrix(0,nrow=ncovs,ncol=ncovs)
for (i in 1:ncovs)
{
   whatread=readBin(z,raw(),16)
   beta.vcv[i,]=readBin(z,numeric(),ncovs,size=8)
}
whatread=readBin(z,raw(),16)
real=list()
real$estimate=readBin(z,numeric(),nlogit,size=8)
whatread=readBin(z,raw(),16)
real$se=readBin(z,numeric(),nlogit,size=8)
whatread=readBin(z,raw(),16)
real$lcl=readBin(z,numeric(),nlogit,size=8)
whatread=readBin(z,raw(),16)
real$ucl=readBin(z,numeric(),nlogit,size=8)
real.vcv=matrix(0,nrow=nlogit,ncol=nlogit)
cont=TRUE
for (i in 1:nlogit)
{
whatread=readBin(z,raw(),16)
whatread=readBin(z,numeric(),nlogit,size=8)
if(length(whatread)!=nlogit)
{
  warning("\n Incomplete read of the binary file\n")
  real.vcv=NULL
  cont=FALSE
  close(z)
  break
} else
{
  real.vcv[i,]=whatread
}
}
if(!is.null(derived_labels)  & cont)
{
	derived=vector("list",length(derived_labels))
	names(derived)=derived_labels
	derived.vcv=vector("list",length(derived_labels))
	names(derived.vcv)=derived_labels
	for(label in derived_labels)
	{
		whatread=readBin(z,raw(),16)
		nderiv=readBin(z,integer(),1,size=4)
		if(length(nderiv)!=0)
		{
			derived[[label]]=list()  
			whatread=readBin(z,raw(),16)
			derived[[label]]$estimate=readBin(z,numeric(),nderiv,size=8)
			whatread=readBin(z,raw(),16)
			derived[[label]]$se==readBin(z,numeric(),nderiv,size=8)
			whatread=readBin(z,raw(),16)
			derived[[label]]$lcl=readBin(z,numeric(),nderiv,size=8)
			whatread=readBin(z,raw(),16)
			derived[[label]]$ucl=readBin(z,numeric(),nderiv,size=8)
			derived[[label]]=as.data.frame(derived[[label]])
			derived.vcv[[label]]=matrix(0,nrow=nderiv,ncol=nderiv)
			for (i in 1:nderiv)
			{
				whatread=readBin(z,raw(),16)
				value=readBin(z,numeric(),nderiv,size=8)
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
