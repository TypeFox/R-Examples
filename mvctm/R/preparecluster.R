preparecluster <-
function(cluster){

cluster=as.matrix(cluster)

N=NROW(cluster)      		 # Total number of observations
nlevel=NCOL(cluster)+1     	 # Number of levels

if(nlevel<2 | nlevel>4){stop("Works for 2, 3 or 4-level data only")}	

id=1:N  			

cl1=unique(cluster[,1])	   	# Labels of the level 1 clusters
m1=length(cl1)  		# Number of level 1 clusters
n1=rep(0,m1)			# Will contain the number of observations in each level 1 cluster

cl2=vector("list",m1)		# Will contain the labels of the level 2 clusters within each level 1 cluster
m2=rep(0,m1)			# Will contain the number of level 2 clusters within each level 1 cluster	
n2=vector("list",m1)		# Will contain the number of observations in each level 2 cluster within each level 1 cluster

cl3=vector("list",m1)		# Will contain the labels of the level 3 clusters within each combination of level 1 and level 2 clusters
m3=vector("list",m1)		# Will contain the number of level 3 clusters nested within each combination of level 1 and level 2 clusters
n3=vector("list",m1) 		# Will contain the number of observations in each level 3 cluster nested within each combination of level 1 and level 2 clusters


clind=vector("list",m1)		# Will be the primary output. The indices to recover the multilevel structure in a nested list form

if(nlevel==2)
	{
	for(i in 1:N)
		{
		ind=which(cl1==cluster[i,1])
		clind[[ind]]=c(clind[[ind]],i) 		
		}
	for(i in 1:m1)	
		{
		n1[i]=length(clind[[i]])
		}
	}

if(nlevel==3)
	{
	for(i in 1:m1)	
		{
		clcur=cluster[cluster[,1]==cl1[i],,drop=F]	# Extracts the level 1 observations
		idcur=id[cluster[,1]==cl1[i]]				# indices of the extracted observations
		cl2[[i]]=unique(clcur[,2])
		n1[i]=NROW(clcur)
		m2[i]=length(cl2[[i]])
		clind[[i]]=vector("list",m2[i])
		n2[[i]]=rep(0,m2[i])
		for(j in 1:n1[i])
			{
			ind=which(cl2[[i]]==clcur[j,2])
			clind[[i]][[ind]]=c(clind[[i]][[ind]],idcur[j]) 		
			}
		}	
	for(i in 1:m1)	
		{
		for(j in 1:m2[i])
			{
		n2[[i]][j]=length(clind[[i]][[j]])
			}
		}
	}


if(nlevel==4)
	{
	for(i in 1:m1)	
		{
		clcur=cluster[cluster[,1]==cl1[i],,drop=F]	# Extracts the level 1 observations
		idcur=id[cluster[,1]==cl1[i]]				# indices of the extracted observations
		n1[i]=NROW(clcur)
		cl2[[i]]=unique(clcur[,2])
		m2[i]=length(cl2[[i]])
		clind[[i]]=vector("list",m2[i])
		n2[[i]]=rep(0,m2[i])
		m3[[i]]=rep(0,m2[i])
		n3[[i]]=vector("list",m2[i])
		cl3[[i]]=vector("list",m2[i])
		for(j in 1:m2[i])
			{
			clcur2=clcur[clcur[,2]==cl2[[i]][j],,drop=F]   # Extracts the level 2 observations within the level 1 cluster
			idcur2=idcur[clcur[,2]==cl2[[i]][j]]		# indices of the extracted observations	
			cl3[[i]][[j]]=unique(clcur2[,3])
			n2[[i]][j]=NROW(clcur2)
			m3[[i]][j]=length(cl3[[i]][[j]])
			clind[[i]][[j]]=vector("list",m3[[i]][j])	
			n3[[i]][[j]]=rep(0,m3[[i]][j])
			for(k in 1:n2[[i]][j])
				{	
				ind=which(cl3[[i]][[j]]==clcur2[k,3])
				clind[[i]][[j]][[ind]]=c(clind[[i]][[j]][[ind]],idcur2[k]) 
				}				
			}		
		}
	for(i in 1:m1)	
		{
		for(j in 1:m2[i])
			{
			for(k in 1:m3[[i]][j])
				{
				n3[[i]][[j]][k]=length(clind[[i]][[j]][[k]])
				}
			}
		}
		
	}

if(nlevel==2){return(list(clind,m1,n1))}
else if(nlevel==3){return(list(clind,m1,n1,m2,n2))}
else if(nlevel==4){return(list(clind,m1,n1,m2,n2,m3,n3))}

}
