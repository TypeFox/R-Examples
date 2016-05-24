makepath <-
# first attempt at making a path for conditional tour, for continuous
# variables only
# Have now added possibility for factors
function (Xc, ncentroids, ninterp = 4)
{
    if (any(apply(Xc, 2L, is.factor))){
        if (!requireNamespace("cluster", quietly = TRUE))
            stop("requires package 'cluster'") 
        d <- cluster::daisy(Xc)
        clustering <- cluster::pam(d, k = ncentroids)  
        centers <- clustering$medoids 
        path <- centers        
    } else {
        if (!requireNamespace("TSP", quietly = TRUE))
            stop("requires package 'TSP'")
        means <- colMeans(Xc)
        sds <- apply(Xc, 2L, sd)
        Xc <- scale(Xc)[, ]
        clustering <- kmeans(Xc, centers = ncentroids)
        centers <- clustering$centers
        o <- TSP::TSP(dist(centers))
        orderindex <- TSP::solve_TSP(o)
        centers <- centers[orderindex, , drop = FALSE]
        interp <- function(x, n = ninterp)
        {
            out <- vector()
            for (i in 1:(length(x) - 1L)){
                out <- c(out, seq(x[i], x[i + 1L], length.out = n + 1L)[-(n + 
                    1L)])
            }
            out
        }
        centers <- as.data.frame(t(apply(t(apply(centers, 1L, `*`, sds)), 1L, 
            `+`, means)))
        path <- as.data.frame(apply(centers, 2L, interp))
    }
    list(centers = centers, path = path)
}