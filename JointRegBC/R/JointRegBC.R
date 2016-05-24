JointRegBC<-function(ini=NA,X,y,z,p,q,...) UseMethod("JointRegBC")
class(JointRegBC)<-"JointRegBC"
JointRegBC.default <- function(ini=NA,X,y,z,p,q, ...)
{
      options(warn = -1)
 f <- function(ini, X, y, z, p, q) {
        X = cbind(1, X)
        y <- as.vector(y)
        z <- as.vector(z)
        ini <- as.vector(ini)
        X <- as.matrix(X)
        n = nrow(X)
        muz = muy = muygivenzx = q2 = q1 = l1 = l2 = l3 = muygivenzx = as.vector(0)
        sez <- ini[p + q + 2]
        seygivenzx <- (1 - (ini[p + q + 1])^2)
mz=matrix(0,n,p)
my=matrix(0,n,q)

for(i in 1:n){
for(j in 1:p){
mz[i,j]=ini[1:p][[j]]*X[i, ][[j]]
}}
for(i in 1:n){
for(k in 1:q){
my[i,k]=ini[(p + 1):(p + q)][[k]]*X[i, -1][[k]]
}}
        for (i in 1:n) {
            muz[i] <- sum(mz[i,])
            muy[i] <- sum(my[i,])
            muygivenzx[i] <- muy[i] + (ini[p + q + 1] * (z[i] - 
                muz[i]))/sez
            q1[i] <- ( - muygivenzx[i])/sqrt(seygivenzx)
            l1[i] <- log(pnorm(q1[i])) + log(dnorm(z[i], muz[i], 
                sez))
            l2[i] <- log(1 - pnorm(q1[i])) + log(dnorm(z[i], 
                muz[i], sez))
        }
        data0 <- cbind(y, l1)
        data1 <- cbind(y, l2)
        data0[data0[, 1] == 1, 2] <- 0
        data1[data1[, 1] == 0, 2] <- 0
        t0 <- sum(data0[, 2])
        t1 <- sum(data1[, 2])
        t <- c(t0, t1)
        Tfinal <- sum(t)
        return(-Tfinal)
    }
    n = nlminb(ini, f, X = X, y = y, z = z, p = p, q = q, lower = c(rep(-Inf, 
        p+q), -0.999, 0), upper = c(rep(Inf, 
        p+q), 0.999, Inf), hessian = T)
    h = fdHess(n$par, f, z = z, y = y, X, p, q)
    h1 = h$Hessian
    ih = ginv(h1)
    se = sqrt(abs(diag(ih)))
    n$Hessian<-h1
    n$p<-p
    n$q<-q
    n$se<-as.vector(se)
    n$call <- match.call()
    class(n) <- "JointRegBC"
object = n
    Co.Re <- data.frame(Parameter = object$par[1:p], S.E = object$se[1:p], 
        `Confidence Interval` = paste("(", round(object$par[1:p] - 
            2 * object$se[1:p], 3), ",", round(object$par[1:p] + 
            2 * object$se[1:p], 3), ")", sep = ""))
    Binary.Re <- data.frame(Parameter = object$par[(p + 1):(p + q)], 
        S.E = object$se[(p + 1):(p + q)], `Confidence Interval` = paste("(", 
            round(object$par[(p + 1):(p + q)] - 2 * object$se[(p + 
                1):(p + q)], 3), ",", round(object$par[(p + 1):(p + 
                q)] + 2 * object$se[(p + 1):(p + q)], 3), ")", 
            sep = ""))
        Cor <- data.frame(Parameter = object$par[p + q + 1], S.E = object$se[p + 
        q + 1], `Confidence Interval` = paste("(", round(object$par[p + 
        q + 1] - 2 * object$se[p + q + 1], 3), ",", round(object$par[p + 
        q + 1] + 2 * object$se[p + q + 1], 3), ")", sep = ""))
    Var <- data.frame(Parameter = object$par[p + q + 2], S.E = object$se[p + 
        q + 2], `Confidence Interval` = paste("(", round(object$par[p + 
        q + 2] - 2 * object$se[p + q + 2], 3), ",", round(object$par[p + 
        q + 2] + 2 * object$se[p + q + 2], 3), ")", sep = ""))
    res <- list(call = object$call, `Continuos Response` = Co.Re, 
        `Variance Of Countinous Response` = Var, `Binary Response` = Binary.Re, 
         Correlation = Cor)
    res$Hessian <- h1
    res$convergence <- n$convergence
res$objective<-  n$objective
    res$call <- match.call()
    class(res) <- "JointRegBC"
    res
}