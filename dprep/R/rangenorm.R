rangenorm <-
function (data, method = c("znorm", "mmnorm", "dscale", "signorm", 
    "softnorm"), superv = TRUE) 
{
    data = as.matrix(data)
    if (!superv) {
        data = cbind(data, rep(1, dim(data)[1]))
    }
    if (!(method %in% c("znorm", "mmnorm", "dscale", "signorm", 
        "softnorm"))) {
        cat("This normalization method  is not supported by this function.\n")
        return(method)
    }
    if (method == "znorm") {
        datanorm = znorm(data)
    }
    if (method == "mmnorm") {
        datanorm = mmnorm(data)
    }
    if (method == "signorm") {
        datanorm = signorm(data)
    }
    if (method == "dscale") {
        datanorm = decscale(data)
    }
    if (method == "softnorm") {
        datanorm = softmaxnorm(data)
    }
    if (!superv) 
        datanorm = datanorm[, -dim(data)[2]]
    return(datanorm)
}
