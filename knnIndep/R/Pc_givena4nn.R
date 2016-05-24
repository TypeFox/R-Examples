Pc_givena4nn <- 
function(i,c,a,k1,k2,N){
	res = rep(-1,length(c))
	if(any(a!=c & k2==1)){
		res[a!=c & k2==1] = Pc_givena(i,c[a!=c & k2==1],a[a!=c & k2==1],N) #k1==0 & k2==1
	}
	if(any(a!=c & k2 != 1)){
			res[a!=c & k2!=1] = (P_cge_aeq(i,c[a!=c & k2!=1],a[a!=c & k2!=1],k2[a!=c & k2!=1],N)-P_cge_aeq(i,c[a!=c & k2!=1]+1,a[a!=c & k2!=1],k2[a!=c & k2!=1],N))/mapply(function(u,v){sum(sapply(u:4,function(l){parameters(l,i-u,v,N)}))},k2[a!=c & k2!=1],a[a!=c & k2!=1])# k1==0 & k2>1
	}
	if(any(a==c)){
		res[a==c] = mapply(function(u,v){sum(sapply(min(4,u+1):4,function(l){parameters(l,i-u,v,N)})) / sum(sapply(u:4,function(l){parameters(l,i-u,v,N)}))},k1[a==c],c[a==c])# k1 > 1
	}
	return(res)
}

