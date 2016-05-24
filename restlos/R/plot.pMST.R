plot.pMST<-function (x, ...) 
{
    # require(rgl)
    drin = x$sample
    data = x$data
    x6 = x$x6
    x6 <- x6[-1, ]
    x7 <- x6[which(duplicated(x6[, 1]) == F), ]
    x7 <- x7[which(x7[, 3] >= floor(dim(data)[1]/2) + 1), ]

    N = length(drin)
    n <- dim(data)[1]
    f <- 1:n
    if (dim(data)[2] == 3) {
        plot3d(data, size = 4)
        points3d(data[drin, ], col = "red", size = 8)
    }
    if (dim(data)[2] == 2) {
        plot(data, xlab = "x", ylab = "y")
        points(data[drin, ], col = "red")
        dev.new()
    }
    if (!is.null(nrow(x7))) {
        par(mfrow = c(2, 2))
    }
    else {
        par(mfrow = c(1, 2))
    }
    plot(mahalanobis(data, colMeans(data[drin, ]), cov(data[drin, 
        ])), xlab = "observation", ylab = "mahalanobis distance", 
        main = "Mahalanobis distances (in-)to subsample")
    points(f[-drin], mahalanobis(data, colMeans(data[drin, ]), 
        cov(data[drin, ]))[-drin], col = "red")
    distl <- as.matrix(dist(data, upper = T, diag = T))
    distl[, -drin] <- rep(max(distl), n)
    diag(distl) <- rep(max(distl), n)
    minIn <- apply(distl[drin, ], 1, min)
    minOut <- apply(distl[-drin, ], 1, min)
    plot(1:n, rep(-1, n), ylim = c(0, max(max(minIn, minOut))), 
        xlab = "observation", ylab = "euclidean distance", main = "Minimal distances (in-)to subsample")
    points(f[drin], minIn)
    points(f[-drin], minOut, col = "red")
	
	if (!is.null(nrow(x7))) {
        plot(x7[, 2], x7[, 1], xlab = "edge length", ylab = "length of largest sub-tree (of MST)", 
            main = "Length-connection-plot (LC plot)")
        lines(x7[, 2], x7[, 1])
        x7 <- cbind(x7, (x7[, 1] - mean(x7[, 1]))/sd(x7[, 1]), 
            (x7[, 2] - mean(x7[, 2]))/sd(x7[, 2]))
        x8 <- as.vector(as.matrix(dist(x7[, 1:2])))
        x9 <- x8[which(x8 == 0) + 1]
        n <- length(x9)
        dx7 <- diff(x7)
        dx7 <- cbind(dx7, dx7[, 4]/dx7[, 5], sqrt(dx7[, 4]^2 + 
            dx7[, 5]^2))
        plot(x7[-1, 3], dx7[, 7]/dx7[, 6], xlab = "obs. in relevant subsample", 
            ylab = "edge-length/ascent ratio (LC plot)", main = "Ascent-length criterion (AL plot)")
        lines(x7[-1, 3], dx7[, 7]/dx7[, 6])
    }
}
