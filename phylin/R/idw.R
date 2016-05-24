idw <-
function(values, coords, grid, method="Shepard", p=2, R=2, N=15) {

    d.real <- real.dist(coords, grid)
    dimensions <- dim(d.real)

    methods <- c("Shepard", "Modified", "Neighbours")
    method <- agrep(method, methods)

    if (method == 1) {
        w <- 1/d.real**p
    } else if (method == 2) {
        w <- ((R-d.real) / (R*d.real))**p
    } else if (method == 3) {
        calcneighbours <- function(x, N) {
            x[order(x)][N:length(x)] <- Inf
            return(x)
        }
        newdist <- t(apply(d.real, 1, calcneighbours, N))
        w <- 1/newdist**p
    }

    # To allow the idw to act on points with same coordinate, rows are checked
    # for infinite weights. When found, points with Inf are 1 and all others 
    # have 0 weight
    for (i in 1:nrow(w)) {
        if (sum(is.infinite(w[i,])) > 0){
            w[i,!is.infinite(w[i,])] <- 0
            w[i,is.infinite(w[i,])] <- 1
        }
    }

    # Interpolation
    w.sum <- apply(w, 1, sum, na.rm=TRUE)
    wx <- w %*% diag(values)
    ux <- apply(wx/w.sum, 1, sum, na.rm=TRUE)

    data.frame(Z=ux)
}
