permute3levelt1 <-
function(cluster,m1,m2){

perm=sample(unlist(cluster,recursive=FALSE))
clusterperm=vector("list",m1)
n2perm=vector("list",m1)

t=1
for(i in 1:m1)
	{
	clusterperm[[i]]=vector("list",m2[i])
	n2perm[[i]]=rep(0,m2[i])
	for(j in 1:m2[i])
		{
		clusterperm[[i]][[j]]=perm[[t]]
		n2perm[[i]][j]=length(perm[[t]])
		t=t+1
		}
	}

list(clusterperm,n2perm)
}
