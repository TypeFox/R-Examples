nonnegreg <-
function(C, S , D)
{
	C <- as.matrix(C)
	S <- as.matrix(S)
	D <- as.matrix(D)
	
	D[is.na(D)] <- 0; C[is.na(C)] <- 0;	S[is.na(S)] <- 0
	
	lin.gain <- max(D)
	D <- D/max(D)
	C <- normalize(C)
	
	C.maxs.pos <- as.vector(apply(C,2,which.max))
	if(ncol(S)!=1) S.max <- t(as.matrix(D[C.maxs.pos,]))
	if(ncol(S)==1) S.max <- as.matrix(D[C.maxs.pos,])

    Sd <- dim(S)
	Sm <- t(apply(as.matrix(1:ncol(D)),1,function(i){
		cofs <- try(nnls(C,D[,i]), silent=T)
		if(class(cofs)!="try-error") {
			cofs <- coefficients(cofs)
			#cofs[cofs<0] <- 0
			}else{
				cofs <- rep(0,ncol(S))
			}
			cofs
		}))
	if(!all(dim(Sm)==Sd)) dim(Sm) <- Sd
	S <- Sm
	S <- apply(as.matrix(1:ncol(S)),1,function(x){apply(matrix(c(S[,x],S.max[,x]),ncol=2),1,min)})
	C <- sweep(C,2,apply(S,2,max),"*")*lin.gain
	S <- normalize(S)
	return(list(C=C,S=S))
}
