mvctm <-

function(fixed,cluster,data,leveltested,method="ls",npermut=1000,weight="observation",affequiv=TRUE)
{

# library("Formula")
# library("nlme")

if(method!="ls" & method!="mixed" & method!="sign" & method!="rank"){stop("Argument method must be one of: ls, mixed, rank, sign")}
if(weight!="observation" & weight!="pair" & weight!="cluster"){stop("Argument weight must be one of: observation, pair, cluster")}
if(!is.data.frame(data)){stop("Argument data must be a data frame")}
if(!is.vector(cluster)){stop("Argument cluster must be a vector")}
if(!(npermut>=1)){stop("Argument npermut must be >= 1")}

clustername=cluster
N=NROW(data) 			# Total number of observations
formu=Formula(fixed)
y=as.matrix(model.part(formu,model.frame(formu,data),lhs=1))
if(NROW(y)!=N){stop("The data cannot contain missing values")}
p=NCOL(y)			# Dimension of the responses
if(method=="mixed" & p>1){stop("Method=\"mixed\" is only available for a univariate response. You could try method=\"ls\".")}

nlevel=length(cluster)+1	# Number of levels

stand="reg"
cova=getCovariateFormula(formu)
if(cova=="~0" | cova=="~-1"){stand="none"}
if(cova=="~1" & method!="mixed"){stand="location"}

cluster=as.matrix(data[,cluster])
if(sum(is.na(cluster))>0){stop("The cluster memberships can not contain missing values")}
if(nlevel!=2 & nlevel!=3 & nlevel!=4){stop("Works for 2, 3 or 4-level data only")}

if(is.numeric(leveltested)==FALSE){stop("You must specify the level to be tested")}
if(length(leveltested)>1){stop("Only one value of leveltested can be given")}
if(nlevel==2){if(leveltested!=1)
		stop("The level to be tested must be 1 for 2-level data")}
else if(nlevel==3)	{if(leveltested!=1 & leveltested!=2)	
		stop("The level to be tested must be 1 or 2 for 3-level data")}
else if(nlevel==4)	{if(leveltested!=1 & leveltested!=2 & leveltested!=3)	
		stop("The level to be tested must be 1, 2 or 3 for 4-level data")}	

# Standardization (centering) of the data prior to performing the permutation test

if(stand=="none")
	{
	yst=y
	smat=diag(p)
	}
else if(stand!="none")
	{
	st=ystand(fixed,clustername,data,method,y,N,p,nlevel,stand,affequiv)	
	yst=st[[1]]
	smat=st[[2]]
	}


# Put the cluster matrix into a usable form for the other functions
prepc=preparecluster(cluster)

if(nlevel==2)
	{
	cluster=prepc[[1]]
	m1=prepc[[2]]
	n1=prepc[[3]]
	pval=testvc2level(yst,cluster,npermut,N,p,smat,m1,n1,weight)
	}

else if(nlevel==3)
	{
	cluster=prepc[[1]]
	m1=prepc[[2]]
	n1=prepc[[3]]
	m2=prepc[[4]]
	n2=prepc[[5]]
	if(leveltested==2)
		{
		pval=testvc3levelt2(yst,cluster,npermut,N,p,smat,m1,n1,m2,n2,weight)
		}
	else if(leveltested==1)
		{
		pval=testvc3levelt1(yst,cluster,npermut,N,p,smat,m1,n1,m2,n2,weight)
		}
	}

else if(nlevel==4)
	{
	cluster=prepc[[1]]
	m1=prepc[[2]]
	n1=prepc[[3]]
	m2=prepc[[4]]
	n2=prepc[[5]]
	m3=prepc[[6]]
	n3=prepc[[7]]	
	if(leveltested==3)
		{
		pval=testvc4levelt3(yst,cluster,npermut,N,p,smat,m1,n1,m2,n2,m3,n3,weight)
		}
	if(leveltested==2)
		{
		pval=testvc4levelt2(yst,cluster,npermut,N,p,smat,m1,n1,m2,n2,m3,n3,weight)
		}
	if(leveltested==1)
		{
		pval=testvc4levelt1(yst,cluster,npermut,N,p,smat,m1,n1,m2,n2,m3,n3,weight)
		}
	}

message1=NA
message2=NA
if(pval==99)		{
	message1="The maximum eigenvalue of the dependency matrix computed with the original data is not positive. The permutation test was not performed. The reported p-value is 1. You could try another method or model."
	pval=1
			}
if(method=="sign" & p==1)	{
	message2="Using method=\"sign\" is not recommended with a univariate response. The test might be liberal."
				}

out=list(pval,messages=list(message1,message2),match.call())
names(out)=c("pvalue","messages","call")

out

}
