#random assignment vector according to randomized block design

blockRand=function(w,nrand,block){
      #formatting
      if(sum(factor(w)==w)>0) w=data.matrix(w)
      if(ncol(w)==1) w=as.vector(w)
      #random assignments
      ctrl = how(within = Within(type = "free"), blocks = as.factor(block))
      if(is.vector(w)){
		set = shuffleSet(length(w), nset = nrand, control = ctrl)
		lapply(1:nrand,function(i) w[set[i,]])
	}else{
		set = shuffleSet(length(w[,1]), nset = nrand, control = ctrl)
		lapply(1:nrand,function(i) w[set[i,],])
	}
}