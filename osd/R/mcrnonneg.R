mcrnonneg <- function(D,C,S,stop.threshold=.0001,max.iter=100)
{	
	
	D[is.na(D)==T] <- 0
	C[is.na(C)==T] <- 0
	S[is.na(S)==T] <- 0
	
	C <- as.matrix(C)
	S <- as.matrix(S)
	
	
	RD <- 1
	res <- rep(0,nrow(D))

	res <- apply(as.matrix(1:nrow(D)),1,function(i) sum((D[i,] - (C[i,] %*% t(S)))^2))
    oldrss <- sum(res)
    
    Cd <- dim(C)
    Sd <- dim(S)
	

	k <- 0
	while(abs(RD)>stop.threshold)
	{
		if(is.even(k)==T){
			S <- t(apply(as.matrix(1:ncol(D)),1,function(i){nnls(C,D[,i])$x}))
			if(!all(dim(S)==Sd)) dim(S) <- Sd
			S <- normalize(S)
		}else{
			C <- t(apply(as.matrix(1:nrow(D)),1,function(i){nnls(S,D[i,])$x}))	
			if(!all(dim(C)==Cd)) dim(C) <- Cd
		}
	
		res <- apply(as.matrix(1:nrow(D)),1,function(i) sum((D[i,] - (C[i,] %*% t(S)))^2))

        rss <- sum(res)
        RD <- ((oldrss - rss)/oldrss)
        RD[is.na(RD)] <- 0
        oldrss <- rss
        k <- k + 1
        if(k>=max.iter) break
	}

	return(list(resC=C,resS=S))
}
