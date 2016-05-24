"rtreeshape" <-
function (n, tip.number, p=0.3, model="", FUN="") {

	if (n==0) {
		return(NULL)
	}

"rtreeshape2" <-
function(tip.number, FUN){
	
	

		
	if (tip.number<2 | tip.number!=floor(tip.number)) {
		stop("tip.number must be an integer greater than 2")
	}
	
	
	merge=matrix(NA, tip.number-1, 2)
	
	merge[1,1]=tip.number
	for (node in 1:(tip.number-1)) {
	
		prob=0:(merge[node,1])
		for (i in 0:(merge[node,1])) {
			prob[i+1]=FUN((merge[node,1]), prob[i+1])
		}
		i=sample(0:(merge[node,1]), 1, prob=prob)
		
		merge[node,2]=i
		if (i!=1) { 
			merge[node+1,1]=i 
		}
		if ((merge[node,1]-i)!=1) { 
			merge[node+i]=(merge[node,1]-i) 
		}
	}
	#names=as.character(1:tip·number)
	res=treebalance(nodes=merge)
	res=as.treeshape(res)
	res
}

Qyule<-function(n,i) {
	if (i==0 | i==n) {0}
	else {1}
}
	
Qaldous<-function(n,i) {
	if (i==0 | i==n) {0}
	else {1/(i*(n-i))}
}

	if (class(tip.number)=='list') {
		tip.number=unlist(tip.number)
	}
	
	if (length(tip.number)>1) {
		res=list()
		current=1
		for (i in 1:length(tip.number)) {
			tmp=rtreeshape(n, tip.number[i], p, model, FUN)
			#if (n==1) {
			#	res[[current]]=tmp
			#	current=current+1
			#} else {
				for (j in 1:n) {
					res[[current]]=tmp[[j]]
					current=current+1
				}
			#}
		}
		return(res)
	} else {
	
		if (n!=floor(n) | n < 0) {
			stop("n must be a positive integer")
		}
	
		if (tip.number!=floor(tip.number) | tip.number < 0) {
			stop("tip.number must be a positive integer")
		}
	
		if (identical(FUN,"")==TRUE & model=="") {stop("at least one option")}
		if (identical(FUN,"")==FALSE & model!="") {stop("at most one option")}
		
		if (identical(FUN,"")==FALSE) {
			#if (n==1) {
			#	return(rtreeshape2(tip.number,FUN))
			#}
			trees=list()
			for (i in 1:n) {
				tree<-rtreeshape2(tip.number,FUN)
				trees[[i]]<-tree
			}
			return(trees)
		}
		
		if (model=="pda"){
			#if (n==1) {
			#	return(rpda(tip.number))
			#}
			trees<-list()
			for (i in 1:n) {
				trees[[i]]<-rpda(tip.number)
			}
			return(trees)
		}
		
		if (model=="yule"){
			#if (n==1) {
			#	return(rtreeshape2(tip.number, Qyule))
			#}
			trees=list()
			for (i in 1:n) {
				#trees[[i]]<-rtreeshape2(tip.number, Qyule)
				trees[[i]]<-ryule(tip.number)
			}
			return(trees)
		}

		if (model=="aldous"){
			#if (n==1) {
			#	return(rtreeshape2(tip.number, Qaldous))
			#}
			trees=list()
			for (i in 1:n) {
				trees[[i]]<-rtreeshape2(tip.number, Qaldous)
			}
			return(trees)
		}


		if (model=="biased"){
			#if (n==1) {
			#	return(rbiased(tip.number=tip.number, p=p))
			#}
			trees<-list()
			for (i in 1:n) {
				trees[[i]]<-rbiased(tip.number=tip.number, p=p)
			}
			return(trees)
		}
		stop("model incorrect")
	}	
}

