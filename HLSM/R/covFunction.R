##input: X is a dataframe with
##	1st column as a vector of network id and is named 'id'
##	2nd column has names of the node

getReceiverCov = function(X,nnodes,nodenames){
	if(all(dimnames(X)[[2]] != 'id')){stop('network id is not provided. Use name id to name the first column of X')}
	if(all(dimnames(X)[[2]] != 'Node')){stop('Nodes are not provided. Use nameNode for the second column')}

	netlist = unique(X$id)
	XX = list()
	PP = dim(X)[2] - 2
	for(ii in 1:length(netlist)){
		XX[[ii]] = array(0, dim = c(nnodes[ii],nnodes[ii],PP))
		Xsm = X[which(X$id == ii),]
		KK = unique(Xsm$Node)
		if(length(KK) != nnodes[ii])(warning('Missing data replaced by 0'))
		dimnames(XX[[ii]])[[1]] = dimnames(XX[[ii]])[[2]] = nodenames[[ii]]
		x1 = which(dimnames(Xsm)[[2]] == 'id')
		x2 = which(dimnames(Xsm)[[2]] == 'Node')
		for(kk in 1:length(KK)){
			k1 = which(nodenames[[ii]] == KK[kk])
			XX[[ii]][,k1,] = sapply(1:PP,function(x)rep(as.numeric(Xsm[which(Xsm$Node==KK[kk]),-c(x1,x2)][x]),nnodes[ii]))
	}
	}
	return(XX)
}


getSenderCov = 	function(X,nnodes,nodenames){
	if(all(dimnames(X)[[2]] != 'id')){stop('network id is not provided. Use name id to name the first column of X')}
	if(all(dimnames(X)[[2]] != 'Node')){stop('Nodes are not provided. Use nameNode for the second column')}

	netlist = unique(X$id)
	XX = list()
	PP = dim(X)[2] - 2
	for(ii in 1:length(netlist)){
		XX[[ii]] = array(0, dim = c(nnodes[ii],nnodes[ii],PP))
		Xsm = X[which(X$id == ii),]
		KK = unique(Xsm$Node)
		if(length(KK) != nnodes[ii])(warning('Missing data replaced by 0'))
		dimnames(XX[[ii]])[[1]] = dimnames(XX[[ii]])[[2]] = nodenames[[ii]]
		x1 = which(dimnames(Xsm)[[2]] == 'id')
		x2 = which(dimnames(Xsm)[[2]] == 'Node')
		for(kk in 1:length(KK)){
			k1 = which(nodenames[[ii]] == KK[kk])
			Xind = which(Xsm$Node==KK[kk])
			for(ll in 1:PP){
			    XX[[ii]][k1,,ll] = rep(as.numeric(Xsm[Xind,
						-c(x1,x2)][ll]),nnodes[ii]) }
	}
	}
	return(XX)
}

getEdgeCov = function(X,nnodes,nodenames){
	if(all(dimnames(X)[[2]] != 'id')){stop('network id is not provided. Use name id to name the first column of X')}
	if(all(dimnames(X)[[2]] != 'Sender')){stop('Sender nodes are not provided. Use name Sender for the second column')}
	if(all(dimnames(X)[[2]] != 'Receiver')){stop('Receiver nodes are not provided. Use name Receiver for the third column')}
	netlist = unique(X$id)
	XX = list()
	PP = dim(X)[2] - 3
	for(ii in 1:length(netlist)){
		XX[[ii]] = array(0, dim = c(nnodes[ii],nnodes[ii],PP))
		Xsm = X[which(X$id == ii),]
		KK = unique(c(Xsm$Receiver,Xsm$Sender))
		if(length(KK) != nnodes[ii])(warning('Missing data replaced by 0'))
		dimnames(XX[[ii]])[[1]] = dimnames(XX[[ii]])[[2]] = nodenames[[ii]]
		KK1 = Xsm$Sender
		KK2 = Xsm$Receiver
		x1 = which(dimnames(Xsm)[[2]] == 'id')
		x2 = which(dimnames(Xsm)[[2]] == 'Sender')
		x3 = which(dimnames(Xsm)[[2]] == 'Receiver')
		for(kk in 1:nrow(Xsm)){		
			k1 = which(nodenames[[ii]] == KK1[kk])
			k2 = which(nodenames[[ii]] == KK2[kk])
			XX[[ii]][k1,k2,] = as.numeric(Xsm[kk,-c(x1,x2,x3)])
	}
	}
	return(XX)
}









