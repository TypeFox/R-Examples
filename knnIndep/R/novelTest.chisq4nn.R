novelTest.chisq4nn <-
		function(xdata,ydata,maxi=length(xdata)-1){
	rx= rank(xdata,ties.method="random")	
	ry= rank(ydata,ties.method="random")	
	N = length(rx)
	
	paths = sapply(1:maxi,generate.paths,rx,ry,N)
	stats = simplify2array(mclapply(1:maxi,function(i){
						if(i==1){
							fi = P_ceq(i,1:floor(N/2),N)
						}else{
							fi =
                            colMeans(sapply(1:floor(N/2),function(c){Pc_givena4nn(i-1,rep(c,maxi),paths[i-1,],rowSums(apply(paths[(i-1):1,,drop=FALSE],1,`==`,c)),rowSums(apply(paths[(i-1):1,,drop=FALSE],1,`==`,paths[i-1,])),N)}))#-Pc_givena4nn(i,rep(c-1,maxi),paths[i-1,],rowSums(apply(paths[(i-1):1,,drop=FALSE],1,`==`,paths[i,])),rowSums(apply(paths[(i-1):1,,drop=FALSE],1,`==`,paths[(i-1),])),N)}))				
						}
						df = sum(fi > 0) -1
						ei = table(factor(paths[i,],levels=1:floor(N/2)))/N
						stat = sum(((ei-fi)^2/fi)[fi>0],na.rm=T)
						return(c(stat,df))
					}),higher=FALSE)
	stat = sum(stats[1,])
	res = list(statistic=stat,p.value=1-pchisq(stat,df=sum(stats[2,])))	
	class(res) = "htest"
	return(res)
}
