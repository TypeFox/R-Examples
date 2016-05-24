flood <-
function ( data, Nx=10, Ny=10, rlen=2000 )
{

if(is.matrix(data)==F)stop("at least two-dimensional data matrix required")
if(mode(data)!="numeric")stop("numeric data required")
if(nrow(data)<=ncol(data))stop("n > d required")

# require(som)


N <- floor((nrow(data)+ncol(data)+1)/2)

rsom <- list()
rsom$som.results <- som(data, xdim=Nx, ydim=Ny, init="linear", alpha=.8, rlen=rlen)



U1 <- matrix(ncol=5, nrow=Nx*Ny)
U1[1:(Nx*Ny), 1] <- 1:(Nx*Ny)
U1[1:(Nx*Ny), 2] <- 1:(Nx*Ny)-1
U1[1:(Nx*Ny), 3] <- 1:(Nx*Ny)+1
U1[1:(Nx*Ny), 4] <- 1:(Nx*Ny)-Nx
U1[1:(Nx*Ny), 5] <- 1:(Nx*Ny)+Nx
U1[U1[, 3]>(sort(rep(c(1:Ny*Nx), Nx))), 3] <- NA
U1[U1[, 2]<(sort(rep(c(1:Ny*Nx-Nx+1), Nx))), 2] <- NA
U1[U1<=0] <- NA
U1[U1>Nx*Ny] <- NA

rsom$som.neigh <- U1



i <- 0
U2 <- numeric(Nx*Ny)
repeat{i <- i+1
	U2[i] <- mean(as.matrix(dist(na.omit(rsom$som.results$code[U1[i, ], ], na.omit=TRUE)))[-1, 1])
	if(i==(Nx*Ny)) break}

U3 <- (U2-min(U2))/(max(U2)-min(U2))
U4 <- matrix(ncol=3, nrow=(Nx*Ny))
U4[, 3] <- U3
U4[, 1] <- rep(c(1:Nx), Ny)
U4[, 2] <- sort(rep(c(1:Ny), Nx))

rsom$umatrix <- U4



T5 <- matrix(1:(Nx*Ny), ncol=Ny, nrow=Nx)
rsom$winneuron <- diag(T5[(rsom$som.results$visual[, 1]+1), (rsom$som.results$visual[, 2]+1)])



T1 <- order(rsom$umatrix[, 3])

l <- 0
LiB <- list()
LiN <- list()
GeB <- c()
FAFH <- matrix(ncol=3, nrow=Nx*Ny)
FAFH_LiB <- list()
FAFH_drin <- list()



repeat{l <- l+1

T2 <- sapply(LiN, function(x){any(x==T1[l])})

if(any(T2)){

	if(sum(T2)>1){
	T4 <- which(T2)
	LiB[[T4[1]]] <- c(unlist(LiB[T4]), T1[l])
	LiB[T4[-1]] <- 0
	LiN[[T4[1]]] <- unique(c(unlist(LiN[T4])), as.vector(na.omit(rsom$som.neigh[T1[l], 2:5])))
	LiN[T4[-1]] <- 0
	GeB[T4[1]] <- sum(GeB[T4])+rsom$som.results$code.sum[T1[l], 3]
	GeB[T4[-1]] <- NA
	}
	else{
	T3 <- which(T2)
	LiB[[T3]] <- c(LiB[[T3]], T1[l])
	LiN[[T3]] <- unique(c(LiN[[T3]], as.vector(na.omit(rsom$som.neigh[T1[l], 2:5]))))
	GeB[T3] <- GeB[T3]+rsom$som.results$code.sum[T1[l], 3]
	}

}
else{
LiB[[l]] <- c(T1[l])
LiN[[l]] <- as.vector(na.omit(rsom$som.neigh[T1[l], 2:5]))
GeB[l] <- rsom$som.results$code.sum[T1[l], 3]
}


FAFH[l, 2] <- max(GeB, na.rm=T)
FAFH[l, 1] <- max(rsom$umatrix[LiB[[which.max(GeB)]], 3])
FAFH[l, 3] <- length(LiB[[which.max(GeB)]])
FAFH_LiB[[l]] <- LiB[[which.max(GeB)]]
FAFH_drin[[l]] <- which(rsom$winneuron%in%LiB[[which.max(GeB)]])

if(any(na.omit(GeB)>=N)==TRUE&exists(x="lib", where=rsom)==FALSE){
	rsom$lib <- LiB
	rsom$lin <- LiN
	rsom$geb <- GeB
	rsom$l <- l
	}

if(l==Nx*Ny) break
}



rsom$fafh <- FAFH
rsom$fafh.lib <- FAFH_LiB
rsom$fafh.drin <- FAFH_drin
rsom$drin <- which(rsom$winneuron%in%rsom$lib[[which.max(rsom$geb)]])

class(rsom) <- "flood"
return(rsom)

}

