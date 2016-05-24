"matsort" <-
function (...) {
    x = cbind(...)
    if(!is.numeric(x))
        stop("Input should by numeric.")
    res = array(dim = c(nrow(x), 0))
    if (!is.matrix(drop(x)))
        return(x)
    else if (ncol(x) > 30)
        return(t(apply(x, 1, sort)))
    else while (is.matrix(drop(x))) {
        imc = max.col(x)
        x = t(x)
        imx = nrow(x) * (1:ncol(x) - 1) + imc
        xmax = x[imx]
        x = t(matrix(x[-imx], ncol = ncol(x)))
        res = cbind(res, xmax)
    }
    res = cbind(res, x)
    colnames(res) = NULL
    rownames(res) = NULL
    return(res)
}
