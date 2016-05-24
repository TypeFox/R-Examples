permute4levelt2 <-
function(cluster,m1,m2,m3){

clusterperm=vector("list",m1)
n3perm=vector("list",m1)
for(i in 1:m1)
	{
	clusteri=cluster[[i]]
	m1i=m2[i]
	m2i=m3[[i]]
	clusterpermi=permute3levelt1(clusteri,m1i,m2i)
	clusterperm[[i]]=clusterpermi[[1]]
	n3perm[[i]]=clusterpermi[[2]]
	}

list(clusterperm,n3perm)
}
