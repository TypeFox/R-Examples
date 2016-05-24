rdela <-
function ( data, N=floor((nrow(data)+ncol(data)+1)/2), rew=TRUE )
{

if(is.matrix(data)==F)stop("at least two-dimensional data matrix required")
if(mode(data)!="numeric")stop("numeric data required")
if(nrow(data)<=ncol(data))stop("n > d required")
if(N<floor((nrow(data)+ncol(data)+1)/2))stop("N >= (n + d + 1)/2 required")

# require(geometry)

alpha <- .975
outp <- list()
outp$data <- data


rdela_tri <- delaunayn(as.matrix(data), options="Fa Fn", full=T)
outp$tri <- rdela_tri$tri

if(any(apply(rdela_tri$tri, 1, function(x) det(cov(data[x, ])))==0))stop("degenerate data set")

rdela_neigh <- lapply(rdela_tri$neighbours, function(x) x[which(x>0)])
outp$neigh <- rdela_neigh

rdela_area <- apply(rdela_tri$tri, 1, function(x) sqrt(-1/sum(solve(as.matrix(dist(data[x, ]))^2*(-.5)))))
outp$radii <- rdela_area
rdela_center <- apply(rdela_tri$tri, 1, function(x) (rowSums(solve(as.matrix(dist(data[x, ]))^2*(-.5)))/sum(solve(as.matrix(dist(data[x, ]))^2*(-.5))))%*%data[x, ])
outp$center <- rdela_center


k <- 0
T1 <- order(rdela_area)
l <- 0
LiB <- list()
LiN <- list()
GeB <- c()
# test88 <- c()
out <- c()
stoppi <- FALSE
cuti <- max(rdela_area)

repeat{l <- l+1

if(rdela_area[[T1[l]]] >= cuti & stoppi) break

T2 <- sapply(LiN, function(x){any(x==T1[l])})

if(any(T2==TRUE)){

	if(sum(T2)>1){
	T4 <- which(T2)
	LiB[[T4[1]]] <- c(unlist(LiB[T4]), T1[l])
	LiB[T4[-1]] <- 0
	LiN[[T4[1]]] <- unique(c(unlist(LiN[T4]), rdela_neigh[[T1[l]]]))
	LiN[T4[-1]] <- 0
	GeB[T4[1]] <- length(unique(as.vector(rdela_tri$tri[LiB[[T4[1]]], ])))
	GeB[T4[-1]] <- NA
	}
	else{
	T3 <- which(T2)
	LiB[[T3]] <- c(LiB[[T3]], T1[l])
	LiN[[T3]] <- unique(c(LiN[[T3]], rdela_neigh[[T1[l]]]))
	GeB[T3] <- length(unique(as.vector(rdela_tri$tri[LiB[[T3]], ])))
	}

}
else{
LiB[[l]] <- c(T1[l])
LiN[[l]] <- rdela_neigh[[T1[l]]]
GeB[l] <- length(unique(as.vector(rdela_tri$tri[LiB[[l]], ])))
}

if(is.null(out) & any(na.omit(GeB)>=N)){
	out <- unique(as.vector(rdela_tri$tri[LiB[[which.max(GeB)]],]))
	em <- unique(as.vector(rdela_area[LiB[[which.max(GeB)]]]))
	m <- length(LiB[[which.max(GeB)]])
	cuti <- mean(em) + sd(em)*min(10, sqrt(abs((m^2-1)/(m^2*(1-alpha)-m))))
	stoppi <- TRUE
	if(rew==FALSE) break
	}

}

outp$LiB <- LiB
outp$LiN <- LiN
outp$GeB <- GeB

outp$drin <- out
outp$raw.mean <- colMeans(data[out,])
outp$raw.cov <- cov(data[out,])
if(rew==TRUE){
outp$final <- unique(as.vector(rdela_tri$tri[LiB[[which.max(GeB)]], ]))
outp$mean <- colMeans(data[outp$final,])
outp$cov <- cov(data[outp$final,])
}

class(outp) <- "rdela"
return(outp)

}

