pred.int.norm.Modified.CA.on.r.prob <-
function (n, df = n - 1, n.mean = 1, K, delta.over.sigma, r, 
    integrate.args.list = NULL) 
{
    if (!is.null(integrate.args.list)) {
        args.list <- c(list(f = pred.int.norm.Modified.CA.on.r.integrand, 
            lower = 0, upper = 1, n.weird = n, df.weird = df, 
            n.mean = n.mean, K.weird = K, delta.over.sigma = delta.over.sigma, 
            r.weird = r), integrate.args.list)
        ret.val <- do.call("integrate", args.list)$value
    }
    else {
        int.list <- integrate(f = pred.int.norm.Modified.CA.on.r.integrand, 
            lower = 0, upper = 1, n.weird = n, df.weird = df, 
            n.mean = n.mean, K.weird = K, delta.over.sigma = delta.over.sigma, 
            r.weird = r, stop.on.error = FALSE)
        if (int.list$message == "OK") {
            ret.val <- int.list$value
        }
        else {
            rel.tol <- 1e-07
            while (rel.tol >= 0.01) {
                int.list <- integrate(f = pred.int.norm.Modified.CA.on.r.integrand, 
                  lower = 0, upper = 1, n.weird = n, df.weird = df, 
                  n.mean = n.mean, K.weird = K, delta.over.sigma = delta.over.sigma, 
                  r.weird = r, rel.tol = rel.tol, stop.on.error = FALSE)
                if (int.list$message == "OK") {
                  ret.val <- int.list$value
                  break
                }
                rel.tol <- rel.tol/10
            }
            if (rel.tol < 0.01) {
                ret.val <- integrate(f = pred.int.norm.Modified.CA.on.r.integrand, 
                  lower = 0, upper = 1, n.weird = n, df.weird = df, 
                  n.mean = n.mean, K.weird = K, delta.over.sigma = delta.over.sigma, 
                  r.weird = r, rel.tol = 0.01)$value
            }
        }
    }
    ret.val
}
