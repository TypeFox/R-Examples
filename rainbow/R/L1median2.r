L1median2 <- function (X, tol = 1e-06, maxstep = 200, na.rm = TRUE, method = c("hossjercroux",
    "coordinate"))
{
    method <- match.arg(method)
    if (method == "coordinate")
        return(apply(X, 2, median.default, na.rm = na.rm))
    else return(hossjercroux(X, tol = tol, maxstep = maxstep,
        na.rm = na.rm))
}
