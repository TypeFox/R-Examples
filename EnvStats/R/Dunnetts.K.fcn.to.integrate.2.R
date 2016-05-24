Dunnetts.K.fcn.to.integrate.2 <-
function (s, d, n, df, k, rho) 
{
    arg.mat <- cbind.no.warn(s = as.vector(s), d = as.vector(d), 
        n = as.vector(n), df = as.vector(df), k = as.vector(k), 
        rho = as.vector(rho))
    for (i in c("s", "d", "n", "df", "k", "rho")) assign(i, arg.mat[, 
        i])
    N <- length(s)
    ret.val <- numeric(N)
    for (i in 1:N) {
        ret.val[i] <- Dunnetts.K.F2(d[i] * s[i], k[i], rho[i]) * 
            dchi(sqrt(df[i]) * s[i], df[i]) * sqrt(df[i])
    }
    ret.val
}
