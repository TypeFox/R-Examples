#bloma <- function(k, ncr = 12, ar = 1){
#	Mlist = list()
#	N=0
#	M=0
#	for(i in 1:k){
#		n <- max(2,rpois(1,ncr))
#		m <- max(2,rpois(1,ar*ncr))
#		Mlist[[i]] <- smat2(n,m,rpois(1,500))
#		N <- N+n
#		M <- M+m
#	}
#	S <- matrix(0,ncol=M,nrow=N)
#	xind <- 1:M
#	yind <- 1:N
#	
#	for(i in 1:k){
#		MS <- Mlist[[i]]
#		
#		
#		xi <- sample(xind, size = ncol(MS))	
#		yi <- sample(yind, size = nrow(MS))
#		S[yi,xi] <- MS
#		xind <- xind[-which(xind %in% xi)]
#		yind <- yind[-which(yind %in% yi)]		
#	}
#	return(S)
#}



getindex = function(ind, dim){
	nd <- length(dim)
	cd <- c(1,cumprod(dim)[-nd])
	return(sum( (ind-1)*cd )+1)
}






arsim <- function(n, dim, k, noise = 0.00, shuffle = TRUE, v = 0.1, minc = 1, exp.prop=NULL, min.prop = 1/dim/4, noise.type = "s", dimnames=list(LETTERS,1:max(dim))){


nd <- length(dim)

stopifnot(all(dim - k*minc >= 0))
if( round((1-noise)*n) == 0){
 k <- 1
 noise <- 0
 }

	if(k==1){
		exp.prop <- 1
		csizes <- matrix(dim, ncol=nd)
	}else{
#clustersizes for each k and dim
csizes <- sapply(dim, function(z){ 
r <- z - k*minc
if(r > 0){
    rmultinom(1,r,rep(1/k,k))
}else{
    rep(0,k)
}
})
csizes <- csizes+minc
	
if(is.null(exp.prop)){
    exp.prop <- (exp.prop<-apply(csizes,1,prod))/sum(exp.prop)
    p1 <- (p1<-runif(k))/sum(p1)
    exp.prop <- (1-v)*exp.prop + v*p1
}
	}
# distribute n
nc <- rmultinom(1,n*(1-noise),prob=exp.prop)

pvecs <- vector("list",nd)

# empty or noisy matrix
#if(noise > 0){
#p1 <- runif(dim[1])
#p2 <- runif(dim[2])
#p1 <- p1/sum(p1)
#p2 <- p2/sum(p2)

#M <- outer(p1,p2)
#if(nd > 2){
#for(i in 3:nd){
#p3 <- 	runif(dim[i])
#p3 <- p3/sum(p3)

#M <- outer(M,p3)	
#}}
#M <- rmultinom(1,round(n*noise,0),M)

#dim(M) <- dim
#}else{
M <- array(0,dim = dim)
#}

# fill in the clusters
lower <- rep(1,nd)
for( s in 1:k ){
p1 <- runif(csizes[s,1])
p2 <- runif(csizes[s,2])
p1 <- p1/sum(p1)
p2 <- p2/sum(p2)
p1 <- p1*(1-min.prop[1])+min.prop[1]
p2 <- p2*(1-min.prop[2])+min.prop[2]
pvecs[[1]] <- c(pvecs[[1]],p1*nc[s]/n)
pvecs[[2]] <- c(pvecs[[2]],p2*nc[s]/n)

CM <- outer(p1,p2)
if(nd > 2){
for(i in 3:nd){
p3 <- 	runif(csizes[s,i])
p3 <- p3/sum(p3)
p3 <- p3*(1-min.prop[i])+min.prop[i]
pvecs[[i]] <- c(pvecs[[i]],p3*nc[s]/n)

CM <- outer(CM,p3)	
}}
CM <- rmultinom(1,nc[s],CM)
upper <- lower + csizes[s,]
spans <- list()
for(j in 1:nd){
spans[[j]] <- lower[j]:(upper[j]-1)
lower[j] <- upper[j]	
}
spans <- expand.grid(spans)
indices <- apply(spans,1,function(z){
getindex(z,dim)	
})
M[indices] <- CM

}
	
if(noise > 0){
	if(noise.type %in% c("I","i","ind","indep")){
		#p1 <- apply(M,1,sum)
		#p2 <- apply(M,2,sum)
        #p1 <- p1/sum(p1)
		#p2 <- p2/sum(p2)
        
        p1 <- pvecs[[1]]
		p2 <- pvecs[[2]]
	
		M2 <- outer(p1,p2)
		if(nd > 2){
			for(i in 3:nd){
				#p3 <- 	apply(M,i,sum)
				#p3 <- p3/sum(p3)
                p3 <- pvecs[[i]]
				M2 <- outer(M2,p3)	
			}}
		M2 <- rmultinom(1,round(n*noise,0),M2)
	
		dim(M2) <- dim
		M <- M+M2
	}
	if(noise.type %in% c("S","s","shuffle")){
		M2 <- rmultinom(1,round(n*noise),M/sum(M))
		dim(M2) <- dim
		ordlist <- list()
		ordlist[[1]] <- M2
		ordlist[2:(nd+1)] <- lapply(dim, function(z) sample(1:z))
		
		M2 <- do.call("[",ordlist)
		M <- M + M2
		}
	}
	if(!is.null(dimnames)){
		if(!is.null(dimnames[[1]])){
			for(i in 1:nd){
				dimnames(M)[[i]] <- rep(dimnames[[1]][i],dim[i])
			}
		}
		if(!is.null(dimnames[[2]])){
			for(i in 1:nd){
				dimnames(M)[[i]] <- paste(dimnames(M)[[i]],dimnames[[2]][1:dim[i]],sep="")
			}
		}
	}
	
	
if( shuffle ){
ordlist <- list()
ordlist[[1]] <- M
ordlist[2:(nd+1)] <- lapply(dim, function(z) sample(1:z))

M <- do.call("[",ordlist)

}

return(as.table(M))
}
