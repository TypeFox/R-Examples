pMST <- function (data, N = floor((nrow(data) + ncol(data) + 1)/2), 
    lmax = nrow(data) * 100) 
{
	if(is.data.frame(data))
		data=as.matrix(data)
		
    if (!is.matrix(data)) 
        stop("at least two-dimensional data matrix required")

    if (mode(data) != "numeric") 
        stop("numeric data required")

    if (dim(data)[1] <= dim(data)[2]) 
        stop("n > d required")
		
    # require(igraph)
    # require(rgl)

    ddmst <- function(dat) {
        ddat <- as.matrix(dist(dat, upper = T, diag = T))
        
		x.tmp <- graph.adjacency(ddat, weighted=TRUE, mode="undirected")
		mstdat  <- minimum.spanning.tree(x.tmp)
		mstdat  <- matrix(as.numeric(get.edgelist(mstdat, names=F)), ncol=2)
		
        o <- dim(dat)[1] - 1
        em <- numeric(o)
        k <- 0
        ddat <- as.matrix(ddat)
        for (i in 1:o) {
            k <- k + ddat[mstdat[i, 1], mstdat[i, 2]]
            em[i] <- ddat[mstdat[i, 1], mstdat[i, 2]]
        }
        emax <- max(em)
        return(k)
    }

    x1 <- as.matrix(dist(data, upper = T, diag = T))
	x.tmp <- graph.adjacency(x1, weighted=TRUE, mode="undirected")
	x2 <- minimum.spanning.tree(x.tmp)

	x2 <- matrix(as.numeric(get.edgelist(x2, names=F)), ncol=2)

    x1 <- as.matrix(x1)
    U1 <- diag(x1[x2[, 1], x2[, 2]])
    T1 <- order(U1)
    l <- 0
    LiB <- list(c())
    GeB <- c()
    LeB <- c()
    x6 <- matrix(c(0, 0, 0), ncol = 3)
    repeat {
        l <- l + 1
        T2 <- sapply(LiB, function(x) {
            any(x == x2[T1[l], 1] | x == x2[T1[l], 2])
        })
        x70 <- 0
        if (any(T2 == TRUE)) {
            if (sum(T2 == TRUE) > 1) {
                T4 <- which(T2 == TRUE)
                maxi <- which.max(sapply(LiB[T4], length))
                x3 <- colMeans(data[unlist(LiB[T4[maxi]]), ])
                LiB[[T4[1]]] <- unique(c(unlist(LiB[T4]), x2[T1[l], 
                  1], x2[T1[l], 2]))
                x4 <- mean(sapply(LiB[T4[-maxi]], function(x) {
                  sqrt(sum((x3 - colMeans(data[x, ]))^2))
                }))
                LiB[T4[-1]] <- 0
                GeB[T4[1]] <- length(LiB[[T4[1]]])
                GeB[T4[-1]] <- NA
                LeB[T4[1]] <- sum(LeB[T4]) + U1[T1[l]]
                LeB[T4[-1]] <- 0
                x5 <- x4
                x7 <- LeB[which.max(sapply(LiB, length))]
            }
            else {
                T3 <- which(T2 == TRUE)
                x3 <- colMeans(data[LiB[[T3]], ])
                LiB[[T3]] <- unique(c(LiB[[T3]], x2[T1[l], 1], 
                  x2[T1[l], 2]))
                x4 <- colMeans(data[x2[T1[l], ], ])
                GeB[T3] <- length(LiB[[T3]])
                x5 <- sqrt(sum((x3 - x4)^2))
                LeB[T3] <- sum(LeB[T3]) + U1[T1[l]]
                x7 <- LeB[which.max(sapply(LiB, length))]
            }
            x6 <- rbind(x6, c(x7, U1[T1[l]], max(sapply(LiB, 
                length))))
        }
        else {
            LiB[[l]] <- c(x2[T1[l], 1], x2[T1[l], 2])
            GeB[l] <- length(LiB[[l]])
            LeB[l] <- U1[T1[l]]
        }
        if (any(na.omit(GeB) >= N) == TRUE || l == lmax) 
            break
    }
    drin <- LiB[[which.max(sapply(LiB, length))]]
    pMST <- list(loc = colMeans(data[drin, ]), cov = cov(data[drin, 
        ]), sample = sort(drin), data = data, x6 = x6)
    class(pMST) <- "pMST"
    return(pMST)
}
