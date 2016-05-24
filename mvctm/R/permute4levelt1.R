permute4levelt1 <-
function(cluster,m1,m2,m3){

perm=sample(unlist(cluster,recursive=FALSE))

clusterperm=vector("list",m1)
n2perm=vector("list",m1)
n3perm=vector("list",m1)
m3perm=vector("list",m1)

t=1
for(i in 1:m1)
	{
	clusterperm[[i]]=vector("list",m2[i])
	n2perm[[i]]=rep(0,m2[i])
	m3perm[[i]]=rep(0,m2[i])
	n3perm[[i]]=vector("list",m2[i])
	for(j in 1:m2[i])
		{
		clusterperm[[i]][[j]]=perm[[t]]
		n2perm[[i]][j]=length(unlist(perm[[t]]))
		m3perm[[i]][j]=length(perm[[t]])
		n3perm[[i]][[j]]=rep(0,m3perm[[i]][j])
		t=t+1
		for(k in 1:m3perm[[i]][j])
			{
			n3perm[[i]][[j]][k]=length(clusterperm[[i]][[j]][[k]])
			}
		}
	}


list(clusterperm,n2perm,m3perm,n3perm)
}
