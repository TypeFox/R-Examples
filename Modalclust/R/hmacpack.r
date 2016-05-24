hmac=function(dat,Sigmas,G=NULL,member=NULL)
{    
	
       dat=as.matrix(dat)
	n=nrow(dat);m=ncol(dat);
	#s_hat=sd(dat);
 	s_hat=apply(dat,2,sd)
       if(is.null(G)) G=dat
       if(is.null(member)) member=seq(n)
        
        
       n.cluster=c();modes=c(); member.n=members=c();
	for (j in 1:length(Sigmas))
	{ if(length(Sigmas)>1) cat("level ",j,"...") 
          else  cat("..........")
          if(j==5) cat("\n")
		sigma=Sigmas[j]; #select j-th bandwidth
		g=1;clust=c();f=c();M=c();density=c()
		if (!is.matrix(G)) G=as.matrix(G)
		#cat("At smoothing level", sigma," the number of clusters are ")
		for (i in 1:nrow(G))
		{                      

                	x0=G[i,,drop=FALSE];x0_old=1;d=1;
                	while (d>10^(-5))
                	{   
			f=mydmvnorm(dat,x0,sigma^2)
                     #f=dmvnorm(dat,x0,diag(sigma^2,m))### still work but slow
                 	p=f/sum(f) #update p
                	x0=p%*%dat; #update x
                    	if (sqrt(sum(x0^2))!=0 ) {d=sqrt(sum((x0_old-x0)^2))/sqrt(sum(x0^2))} 	
			else d=10; #stopping rule
                    	x0_old=x0;
                	}
			
			if (length(dim(M))!=0)
			{
			if(m==1) temp=c(1:dim(M)[1])[abs(M[,1]-x0[1])<s_hat[1]*0.001]
           		else 	  temp=c(1:dim(M)[1])[abs(M[,1]-x0[1])<s_hat[1]*0.001 & abs(M[,2]-x0[2])<s_hat[2]*0.001];
              	if (length(temp)==0) {clust[i]=g;g=g+1} 
			else {clust[i]=clust[temp[1]]}
          		}
 
		else {clust[i]=g;g=g+1}
		M=rbind(M,x0);
		M=as.matrix(M)
             	member.n[member==i]=clust[i]           
		}           	
		G=M[match(c(1:max(clust)),clust),]

            	#save # of clusters by each bandwidth
            	n.cluster[j]=max(clust);
            	#save the variables
            	modes[[j]]=G
            	members[[j]]=member.n
            	#update membership 
            	member=member.n;           
    	}

       output=list() 
	output[["data"]]=dat
    	output[["n.cluster"]]=n.cluster
    	output[["level"]]=as.vector(cumsum(!duplicated(n.cluster)))
	output[["mode"]]=as.vector(modes[c(1:length(Sigmas))[!duplicated(n.cluster)]])
    	output[["membership"]]=as.vector(unique(members))
	class(output)="hmac"
	return(output)

}



phmac=function(dat,length=10,npart=1,parallel=TRUE,sigmaselect=NULL,G=NULL)
{
if(.Platform$OS.type=="windows") {
  parallel=FALSE
  cat ("Multicore processing not implemented in Windows version\n")
}
### Sigma  
  	#scale.dat=sd(dat)
  	scale.dat=apply(dat,2,sd)
	#Data=scale(dat,scale=sd(dat),center=FALSE)
	Data=scale(dat,scale=scale.dat,center=FALSE)
	
    	
    	if (!is.matrix(Data)) Data=as.matrix(Data)
       n=dim(Data)[1];
       m=dim(Data)[2];

	if(m==1) { Data=sort(Data) 
	if(npart!=1) npart=1
	warning("partitioning one dimensional data may cause the clustering to procedure work inproperly.  
	This program has forced the data being sorted and using nonpartitioning option")
	}
      	Data=as.matrix(Data)

       if(is.null(sigmaselect)) Sigmas=khat.inv(p=m,len=length) ### NO sigma given
       else Sigmas=sigmaselect;
	
	if(is.null(G)) G=Data
       n.cluster=c();member=seq(n);member.n=c();modes=c();members=c()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  %%%%%%%  Start here
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Call hmac to do partition on part of the data
        partition=sample(seq(1:npart),nrow(Data),replace=TRUE) #### Can be optimized
	 #partition=rep(1:npart,length=n)
        data.split=split(data.frame(Data),partition)
        if(npart>1) cat("Partioning Data\n")
        else cat("Performing initial Modal clustering")
        if(parallel){
        
          if(require(parallel)){
            if(npart>1) cat("Using parallel computing for performing initial Modal clustering \n")
            split.hmac <- mclapply(1:npart, function(x) hmac(data.split[[x]],Sigmas=Sigmas[1]), mc.cores=npart)
                      }
          else{
            cat("Use library parallel to run hmac in several processors ")  
            split.hmac=lapply(data.split,hmac,Sigmas=Sigmas[1])
          }
        }
        else split.hmac=lapply(data.split,hmac,Sigmas=Sigmas[1])
        
 
        	G=NULL;
      		member.n=rep(0,n);


        if(npart>1)  cat("\nBuilding membership from initial partitions")

        start=0
        for( i in 1:npart){
          if(m >1) G=rbind(G,split.hmac[[i]]$mode[[1]])  #for1d
          else G=c(G,split.hmac[[i]]$mode[[1]])
          member.n[partition==i]=start+split.hmac[[i]]$membership[[1]]
          start=start+split.hmac[[i]]$n.cluster
          
        }
	member=member.n;
        cat("\nBuilding hierarchical Modal clusters at \n") 
      	output= hmac(Data,Sigmas,G,member.n) 
      
	for(j in 1:length(unique(output$level) )){
	output$mode[[j]]=scale(matrix(output$mode[[j]],ncol=m),scale=1/scale.dat,center=FALSE)}
	output$data=dat
	output$sigmas=Sigmas
       cat("\n\n")
	return(output)
 
}


modalclust <- phmac


### This program will be used to obtain the correct Sigma's fpr model clustering
sdofnorm<-function(h,p){
  num= (h^(-p)-(h^2+2)^(-p/2))^2
   din=  (h^(-p)*(h^2+4)^(-p/2))-2*((h^2+1)^(-p/2)*(h^2+3)^(-p/2))+(h^2+2)^(-p)
    dof=        num/din
approx=(sqrt((h^2+4)/h^2))^p
  return(cbind(dof,approx))
}



##########################################
khat<-function(dof,p)
  return(p*(dof^(1/p)-1)/2)

##########################################
khat.inv<- function(p,len=10){	
	mygrid=seq(1+.4*p,6+.4*p,length=len)
	khat.vec=(seq(0,20,by=.01))
  
	hvec=2/(1+2*khat.vec/p)
	for(i in 1:length(khat.vec))
	khat.vec[i]=khat(sdofnorm(hvec[i],p)[,1],p)

	mygrid.h=mygrid
	for(i in 1:length(mygrid)){
		mygrid.h[i]=mean(hvec[abs(khat.vec-mygrid[i])<1e-2])
	}
	return(rev(mygrid.h))
}


mydmvnorm <- function (x, mean, sigmasq){ 
# sigmasq is common variance
	x=as.matrix(x)
    	distval <- mahalanobis(x, center = mean, cov = diag(1/sigmasq,ncol(x)),inverted=TRUE)
    	logretval <- -(ncol(x)*(log(2 * pi) + log(sigmasq)) + distval)/2
    	return(exp(logretval))
}
