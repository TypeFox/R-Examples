bic.generation <-
function(gene.data,affy.loc,set.data)
{
	#====  UnitTests for the function bic.generation
	num.gene.data=sum(apply(is.finite(gene.data),1,sum))
	if(num.gene.data!=dim(gene.data)[1]*dim(gene.data)[2])
		stop("Some entries in the gene.data matrix ar nor numerical")	
	num.set.data=sum(apply(is.finite(set.data),1,sum))
	if(num.set.data!=dim(set.data)[1]*dim(set.data)[2])
		stop("Some entries in the gene.set data matrix ar nor numerical")
	num.afyy.loc=intersect(affy.loc,c(1:dim(gene.data)[1]))
	if(length(num.afyy.loc)!=length(affy.loc))
		stop("Some entries of the affy identifiers do not correspond to the gene.set data matrix")	
		
	#========================================================================	
	bic.gene.0=matrix(0,length(affy.loc),ncol=1)
	bic.gene.1=matrix(0,length(affy.loc),length(affy.loc))
	aa=10
	percent.10=floor(0.1*length(affy.loc))
	for(i in 1:(length(affy.loc)-1))
	{
		modelo=lm(gene.data[affy.loc[i],]~1)
		bic.gene.0[i]=BIC(modelo)
		aux=lapply(((i+1):length(affy.loc)),function(x){
	              modelo=lm(gene.data[affy.loc[i],]~gene.data[affy.loc[x],])
	              BIC(modelo)})    
    	bic.gene.1[i,((i+1):length(affy.loc))]=unlist(aux)    
    	if(i>1)
    	{
    		aux=lapply((1:(i-1)),function(x){
	              modelo=lm(gene.data[affy.loc[i],]~gene.data[affy.loc[x],])
	              BIC(modelo)})
	    	bic.gene.1[i,(1:(i-1))]=unlist(aux)        
    	}
    	if(i%%percent.10==0)
    	{
    		print("BICs in genes completed:")
    		print(paste(aa,"%",sep=""))
    		aa=aa+10
    	}
	}
	modelo=lm(gene.data[affy.loc[length(affy.loc)],]~1)
	bic.gene.0[length(affy.loc)]=BIC(modelo)
	aux=lapply((1:(length(affy.loc)-1)),function(x){
	              modelo=lm(gene.data[affy.loc[length(affy.loc)],]~gene.data[affy.loc[x],])
	              BIC(modelo)})
	bic.gene.1[length(affy.loc),(1:(length(affy.loc)-1))]=unlist(aux)
	
	aa=10
	percent.10=floor(0.1*dim(set.data)[1])
	bic.fc.0=matrix(0,dim(set.data)[1],ncol=1)
	bic.fc.1=matrix(0,dim(set.data)[1],dim(set.data)[1])
	for(i in 1:(dim(set.data)[1]-1))
	{
		modelo=lm(set.data[i,]~1)
		bic.fc.0[i]=BIC(modelo)
		aux=lapply(((i+1):dim(set.data)[1]),function(x){
	              modelo=lm(set.data[i,]~set.data[x,])
	              BIC(modelo)})
		bic.fc.1[i,((i+1):dim(set.data)[1])]=unlist(aux)
		if(i>1)
    	{
    		aux=lapply((1:(i-1)),function(x){
	              modelo=lm(set.data[i,]~set.data[x,])
	              BIC(modelo)})
	     	bic.fc.1[i,(1:(i-1))]=unlist(aux)        
    	}
    	if(i%%percent.10==0)
    	{
    		print("BICs in gene sets completed:")
    		print(paste(aa,"%",sep=""))
    		aa=aa+10
    	}
	}
	modelo=lm(set.data[dim(set.data)[1],]~1)
	bic.fc.0[dim(set.data)[1]]=BIC(modelo)
	aux=lapply((1:(dim(set.data)[1]-1)),function(x){
	              modelo=lm(set.data[dim(set.data)[1],]~set.data[x,])
	              BIC(modelo)})
	bic.fc.1[dim(set.data)[1],(1:(dim(set.data)[1]-1))]=unlist(aux)	
	
	result=list(bic.gene.0,bic.gene.1,bic.fc.0,bic.fc.1)
	nombres=c("BIC.Gene.0.Pred","BIC.Gene.1.Pred","BIC.Set.0.Pred","BIC.Set.1.Pred")
	names(result)=nombres
	return(result)
}
