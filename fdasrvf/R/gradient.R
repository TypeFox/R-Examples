#' Gradient using finite diferences
#'
#' This function takes the gradient of f using finite differences
#'
#' @param f vector with \eqn{N} samples
#' @param binsize scalar of time samples
#' @return g vecotr with \eqn{N} samples which is the gradient of f
#' @keywords srvf alignment
#' @export
#' @examples
#' data("simu_data")
#' out = gradient(simu_data$f[,1],mean(diff(simu_data$time)))
gradient <- function(f,binsize){
    n = nrow(f)
    if (is.null(n)){
        f = as.vector(f)
        n = length(f)
        g = rep(0,n)
        h = binsize*(1:n)
        # Take forward differences on left and right edges
        g[1] = (f[2] - f[1])/(h[2]-h[1])
        g[n] = (f[n] - f[(n-1)])/(h[length(h)]-h[(length(h)-1)])

        # Take centered differences on interior points
        h = h[3:n]-h[1:(n-2)]
        g[2:(n-1)] = (f[3:n]-f[1:(n-2)])/h[1]
    }else {
        f = as.matrix(f)
        p = ncol(f)
        g = matrix(0,n,p)
        h = binsize*(1:n)
        # Take forward differences on left and right edges
        g[1,] = (f[2,] - f[1,])/(h[2]-h[1])
        g[n,] = (f[n,] - f[(n-1),])/(h[length(h)]-h[(length(h)-1)])

        # Take centered differences on interior points
        h = h[3:n]-h[1:(n-2)]
        g[2:(n-1),] = (f[3:n,]-f[1:(n-2),])/h[rep(1,p)]
    }

    return(g)
}
