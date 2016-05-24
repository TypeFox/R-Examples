maxempty <- function(x, y, ax = c(0, 1), ay = c(0, 1)) {

    n <- length(x)
    d <- sort(c(ax, x))
    D <- diff(d)
    m <- which.max(D)
    # check vertical slices
    mgap <- D[m]
    maxr <- mgap * (ay[2] - ay[1])
    maxR <- c(d[m], ay[1], d[m+1], ay[2])
    o <- order(y)
    X <- x[o]; Y <- y[o]
    for (i in 1:n) {
        tl <- ax[1]; tr <- ax[2]
        if (i < n) {
            for (j in (i+1):n) {
                if (X[j] > tl && X[j] < tr) {
                    # check horizontal slices (j == i+1)
                    # and (all) rectangles above (X[i], Y[i])
                    area <- (tr-tl)*(Y[j]-Y[i])
                    if (area > maxr) {
                        maxr <- area
                        maxR <- c(tl, Y[i], tr, Y[j])
                    }
                    if (X[j] > X[i]) tr <- X[j]
                    else             tl <- X[j]
                }
            }
        }
        # check open rectangles above (X[i], Y[i])
        area <- (tr-tl)*(ay[2]-Y[i])
        if (area > maxr) {
            maxr <- area
            maxR <- c(tl, Y[i], tr, ay[2])
        }
        # check open rectangles below (X[i], Y[i])
        ri <- min(ax[2], X[Y < Y[i] & X > X[i]])
        li <- max(ax[1], X[Y < Y[i] & X < X[i]])
        area <- (ri-li)*(Y[i]-ay[1])
        if (area > maxr) {
            maxr <- area
            maxR <- c(li, ay[1], ri, Y[i])
        }
    }
    return(list(area = maxr, rect = maxR))
}
