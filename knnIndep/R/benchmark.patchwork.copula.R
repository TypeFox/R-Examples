benchmark.patchwork.copula <-
function(fun #function to benchmark
		,args #list of arguments for the function fun
		,cvals #,target.mis # mutual information values
		,n=320,nsim=500
		,bins=20
){
	copula.vals = mclapply(1:nsim,function(x){
				lapply(cvals,function(cval){
							dep = generate.patchwork.copula(c=cval,returnmi=TRUE,p=matrix(rbeta(bins*bins,.01,1),ncol=bins),bins=bins,npoints=n,plot=FALSE)
							x2 = runif(n)
							q = do.call(fun,args=c(list(x2,dep$y),args))
							s = do.call(fun,args=c(dep[c("x","y")],args))
							return(c(q,s,mi=dep$mi))
						})
			})
	
	l = lapply(copula.vals,function(runs){lapply(runs,function(vals){list(q=vals[1],s=vals[2])})})
	arr = array(unlist(l,recursive=TRUE),dim=c(500,7,2))
	vals = lapply(1:dim(arr)[1],function(i){lapply(1:dim(arr)[3],function(j){arr[i,,j]})})
	vals = list(lapply(vals,function(lis){names(lis)<-c("null","dep");return(lis)}))
	return(vals)
}
