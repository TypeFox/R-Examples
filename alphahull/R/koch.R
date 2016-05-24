koch <-
function (side = 3, niter = 5) 
{
    npoints = 3 * 4^(niter - 1)
    p <- matrix(nrow = npoints, ncol = 2)
    index <- c(1, npoints/3 + 1, 2 * npoints/3 + 1)
    p[1, ] <- c(-side/2, 0)
    p[index[3], ] <- c(side/2, 0)
    p[index[2], ] <- p[1, ] + rotation(p[index[3], ] - p[1, ], 
        -pi/3)
    if (niter >= 2) {
        for (k in 2:niter) {
            npoints2 = 3 * 4^(niter - k)
            index1 <- index + npoints2/3
            index2 <- index + npoints2
            index3 <- index1 + npoints2/3
            aux <- diff(rbind(p[!is.na(p[, 1]), ], p[1, ]))/3
            p[index1, ] <- p[index, ] + aux
            p[index2, ] <- p[index, ] + 2 * aux
            aux1 <- matrix(c(cos(-pi/3), sin(-pi/3)), nrow = length(index1), 
                ncol = 2, byrow = TRUE)
            aux2 <- matrix(c(-sin(-pi/3), cos(-pi/3)), nrow = length(index1), 
                ncol = 2, byrow = TRUE)
            v <- p[index2, ] - p[index1, ]
            p[index3, ] <- p[index1, ] + cbind(apply(aux1 * v, 
                1, sum), apply(aux2 * v, 1, sum))
            index <- sort(c(index, index1, index2, index3))
        }
    }
    return(vertices = p)
}
