qrnn.rbf <-
function(x, x.basis, sigma)
{
    kern <- matrix(0, nrow=nrow(x), ncol=nrow(x.basis))
    for (k in 1:nrow(x.basis)){
        x.basis.test <- matrix(x.basis[k,], nrow=nrow(x),
                               ncol=ncol(x.basis), byrow=TRUE)
        kern[,k] <- exp(-apply(((x-x.basis.test)^2)/(2*sigma^2), 1, sum))
    }
    kern
}
