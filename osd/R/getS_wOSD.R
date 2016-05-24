getS_wOSD <-
function(C,S,D)
{
	#Cs <<- C
	#Ds <<- D
	#Ss <<- S
	
	C <- as.matrix(C)
	S <- as.matrix(S)
	D <- as.matrix(D)

	D[is.na(D)] <- 0; C[is.na(C)] <- 0;	S[is.na(S)] <- 0
	
	Sd <- dim(S)
	S.reg <- S*0
	C.n <- C*0
	S.reg <- apply(as.matrix(C), 2, function(mod.c) getS.OSD(mod.c,D))
	if(!all(dim(S.reg)==Sd)) dim(S.reg) <- Sd

	# apex.h <- apply(as.matrix(S.reg),2,max)
	# max.frag <- as.vector(apply(S.reg,2,which.max))
	
	# apex.h <- apply(as.matrix(1:length(max.frag)),1,function(i){
		# cofs <- try(nnls(as.matrix(normalize(C[,i])),D[,max.frag[i]]), silent=T)
		# if(class(cofs)!="try-error") {
			# #cofs <- coefficients(cofs)[-1]
			# cofs <- coefficients(cofs)
			# #cofs <- cofs[i]
			# if(is.na(cofs)) cofs <- 0
			# if(cofs<0) cofs <- 0
			# if(cofs>max(D[,max.frag[i]])) cofs <- max(D[,max.frag[i]])
			# }else{
				# cofs <- 0
			# }
			# cofs
		# })
	#C.n <- sweep(normalize(C),2,apex.h,"*")

	S.reg <- normalize(S.reg)
	C <- t(apply(as.matrix(1:nrow(D)),1,function(i){nnls(normalize(S.reg),D[i,])$x}))	
	if(ncol(C)!=ncol(S.reg)) C <- t(C)
	
	if(ncol(C)!=1)
	{
		Cmteo <- apply(D[as.vector(apply(C,2, which.max)),],1,max)
		Cmemp <- apply(C,2,max)
		adj.q <- Cmteo/Cmemp
	}else{
		Cmteo <- max(D[as.vector(apply(C,2, which.max)),])	
		Cmemp <- max(C)
		adj.q <- Cmteo/Cmemp
	}
	adj.q[adj.q>1] <- 1
	C <- sweep(C,2,adj.q,"*")	
	
		
	return(list(C=C, S=S.reg))		
}
