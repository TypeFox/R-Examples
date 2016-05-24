run.tests <-
function(fun, args,types
		, noises #vector or Nxtypes matrix with N noise levels for each relationship type
		, size = 320 #number of datapoints
		,nsim=500 # The number of null datasets we use to estimate our rejection reject regions for an alternative with level 0.05
		,... #optional arguments to function generate.benchmark.data
){
	vals  = vector("list",length(types))
	for(typ in types){
		samples = lapply(1:nsim,function(i){generate.benchmark.data(typ,noises[,which(types == typ)],size,...)})
		vals[[which(types == typ)]] = mclapply(samples,function(xylist){
					val = sapply(1:ncol(xylist$x),function(col){do.call(fun,args=c(list(runif(size),xylist$y[,col]),args))})
					val2 = sapply(1:ncol(xylist$x),function(col){do.call(fun,args=c(list(xylist$x[,col],xylist$y[,col]),args))})
					return(list(null=val,dep=val2))
				})
		
	}
	return(vals)
}
