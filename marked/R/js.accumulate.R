#' Accumulates common capture history values
#' 
#' To speed up compuation, animals with the same capture history and design matrix are
#' accumulated and represented by a frequency. Computes starting values for Phi and p parameters from the
#' list of design matrices and the summarized data list including ch matrix and
#' first and last vectors. If any values are missing (NA) or abs(par)>5, they are
#' set to 0.
#' 
#' 
#' @param x data
#' @param model_data list of design matrices, fixed parameters and time intervals all which can vary by animal
#' @param nocc number of capture occasions
#' @param freq frequency of each capture history before accumulation
#' @param chunk_size size that determines number of pieces of data/design matrix that are handled. Smaller chunk_size requires more
#' time but less memory. 1e7 is default set in cjs. 
#' @return modified model_data list that is accumulated
#' @author Jeff Laake
#' @keywords utility
js.accumulate=function(x,model_data,nocc,freq,chunk_size)
{	
chmat=model_data$imat$chmat
if(sum(model_data$imat$loc)>0)
{
	indices=which(model_data$imat$loc==1)
	chmat[cbind(indices,model_data$imat$last[indices])]=2
	ch=apply(chmat,1,paste,collapse="")
} else
    ch=x$ch
#     Create fixed parameter matrices
Phifixmat=NULL
pfixmat=NULL
pentfixmat=NULL
if(model_data$Phi.fixed[1,1]>0)
{
	Phifixmat=rep(NA,nrow(x)*(nocc-1))
	Phifixmat[(nocc-1)*(model_data$Phi.fixed[,1]-1)+model_data$Phi.fixed[,2]]=model_data$Phi.fixed[,3]      
}
if(model_data$p.fixed[1,1]>0)
{         
	pfixmat=rep(NA,nrow(x)*nocc)
	pfixmat[nocc*(model_data$p.fixed[,1]-1)+model_data$p.fixed[,2]-1]=model_data$p.fixed[,3]      
}
if(model_data$pent.fixed[1,1]>0)
{         
	pentfixmat=rep(NA,nrow(x)*(nocc-1))
	pentfixmat[(nocc-1)*(model_data$pent.fixed[,1]-1)+model_data$pent.fixed[,2]-1]=model_data$pent.fixed[,3]      
}
#     Start accumulating 
#     Construct common splits based on Phi and pent dms and fixed parameter values
nrows=nrow(model_data$Phi.dm)
pieces=floor(max(ncol(model_data$Phi.dm),ncol(model_data$pent.dm))*nrows/chunk_size)+1
chdesign=NULL
iseq=seq(1,nrow(x),floor(nrows/pieces/(nocc-1)))
if(iseq[length(iseq)]!=nrow(x))iseq=c(iseq,nrow(x))
prev_i=1
for(i in iseq[-1])
{
	lower=(prev_i-1)*(nocc-1)+1
	upper=i*(nocc-1)
	xdesign=cbind(rep(ch[prev_i:i],each=nocc-1),as.matrix(model_data$Phi.dm[lower:upper,,drop=FALSE]),as.matrix(model_data$pent.dm[lower:upper,,drop=FALSE]),
			Phifixmat[lower:upper],pentfixmat[lower:upper])
	xdesign=sapply(split(xdesign,rep(1:(i-prev_i+1),each=nocc-1)),paste,collapse="")
	chdesign=c(chdesign,xdesign)
	prev_i=i+1
}
#     If time intervals vary across individuals, then use it as part of accumulation 
#     split
occ.int=as.vector(unlist(apply(model_data$time.intervals,2,unique)))
if(length(occ.int)!=(nocc-1))
	chdesign=paste(chdesign,apply(model_data$time.intervals,1,paste,collapse=""),sep="")
#     Now use p which has different parameter structure	 	  
nrows=nrow(model_data$p.dm)
pieces=floor(ncol(model_data$p.dm)*nrows/chunk_size)+1
chdesignp=NULL
iseq=seq(1,nrow(x),floor(nrows/pieces/(nocc)))
if(iseq[length(iseq)]!=nrow(x))iseq=c(iseq,nrow(x))
prev_i=1
for(i in iseq[-1])
{
	lower=(prev_i-1)*nocc+1
	upper=i*nocc
	xdesign=cbind(rep(ch[prev_i:i],each=nocc),as.matrix(model_data$p.dm[lower:upper,,drop=FALSE]),
			pfixmat[lower:upper])
	xdesign=sapply(split(xdesign,rep(1:(i-prev_i+1),each=nocc)),paste,collapse="")
	chdesignp=c(chdesignp,xdesign)
	prev_i=i+1
}
#     Split data based on paste of all dms, fixed values
chsplit=split(1:nrow(x),paste(chdesign,chdesignp,as.numeric(freq==0),sep=""))
#     Get new set of indices
indices=as.vector(sapply(chsplit,min))
#     Accumulate frequencies
counts=as.vector(sapply(chsplit,function(x)sum(freq[x])))
freq=counts[order(indices)]
#     Get reduced set of dms and fixed parameter matrices
dm.index=rep(1:nrow(x),each=nocc-1)%in%indices
model_data$Phi.dm=model_data$Phi.dm[dm.index,,drop=FALSE]
model_data$pent.dm=model_data$pent.dm[dm.index,,drop=FALSE]
if(!is.null(Phifixmat))
{
	chi=(rep(1:nrow(x),each=nocc-1))[dm.index][!is.na(Phifixmat)[dm.index]]
	occj=(rep(1:(nocc-1),times=nrow(x)))[dm.index][!is.na(Phifixmat)[dm.index]]
	val=Phifixmat[dm.index][!is.na(Phifixmat)[dm.index]]
	model_data$Phi.fixed=cbind( match(chi,indices),occj,val)
} 
if(!is.null(pentfixmat))
{
	chi=(rep(1:nrow(x),each=nocc-1))[dm.index][!is.na(pentfixmat)[dm.index]]
	occj=(rep(1:(nocc-1),times=nrow(x)))[dm.index][!is.na(pentfixmat)[dm.index]]
	val=pentfixmat[dm.index][!is.na(pentfixmat)[dm.index]]
	model_data$pent.fixed=cbind( match(chi,indices),occj,val)
}  
dm.index=rep(1:nrow(x),each=nocc)%in%indices
model_data$p.dm=model_data$p.dm[dm.index,,drop=FALSE]
if(!is.null(pfixmat))
{
	chi=(rep(1:nrow(x),each=nocc))[dm.index][!is.na(pfixmat)[dm.index]]
	occj=(rep(1:nocc,times=nrow(x)))[dm.index][!is.na(pfixmat)[dm.index]]
	val=pfixmat[dm.index][!is.na(pfixmat)[dm.index]]
	model_data$p.fixed=cbind( match(chi,indices),occj,val)
}  
ch=x$ch[sort(indices)]
#     Compute new imat
model_data$imat=process.ch(ch,freq,all=TRUE)
model_data$time.intervals=model_data$time.intervals[sort(indices),]
nx=sum(x$freq)
if(sum(freq)!=nx)
	stop(paste("Error in accumulation. Number of accumulated",sum(freq),"not equal to original number",nx))
cat(" ",nrow(x)," capture histories collapsed into ",length(ch),"\n")
return(model_data)
}

