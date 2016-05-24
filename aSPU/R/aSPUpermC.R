## Sum of powered score (SPU) test (perm, permutation part coded in C)
##
## It gives the p-values of the SPU test and aSPU test based on based on the permutation of residuals.
##
## @param Y phenotype data. It can be disease lables; =0 for controls, =1 for cases.
##     or It can be any quantitative traits. Vector with length n (number of observations)
##
## @param X genotype data; each row for a subject, and each column
##     for an SNP. The value of each element is the # of the copies
##     for an allele. Matrix with dimension n by k (n : number of observation, k : number of genotype data)
##
## @param cov covariates. Matrix with dimension n by p (n :number of observation, p : number of covariates)
##
## @param model Use "gaussian" for quantitative trait (Default)
##    , and Use "binomial" for binary trait.
##
## @param pow power used in SPU test. Vector of g number of power.
##
## @param n.perm number of permutation
##
## @param userank similar to the original code if TRUE, by definition if FALSE
##
## @export
## @return Test Statistics and p-values for SPU tests and aSPU test.
##
## @examples
##
## data(exdat)
## out <- aSPUpermC(exdat$Y, exdat$X, cov = NULL,
##               model = "binomial", pow = c(1:8, Inf), n.perm = 1000)
## out
##
## @seealso \code{\link{aSPU}}


aSPUpermC <- function(Y, X, cov = NULL, model=c("gaussian","binomial"), pow=c(1:8, Inf), n.perm=1000, userank = T){

    model = match.arg(model)

    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)

#### Score vector:
    if (is.null(cov)){
        ## NO nuisance parameters:
        XUs<-Xg <- X
        r<-Y-mean(Y)
        U<-as.vector(t(Xg) %*% r)
   } else {
       tdat1<-data.frame(trait=Y, cov)
       fit1<-glm(trait~.,family=model,data=tdat1)
       pis<-fitted.values(fit1)
       XUs<-matrix(0, nrow=n, ncol=k)
       Xmus = X
       for(i in 1:k){
           tdat2<-data.frame(X1=X[,i], cov)
           fit2<-glm(X1~.,data=tdat2)
           Xmus[,i]<-fitted.values(fit2)
           XUs[, i]<-(X[,i] - Xmus[,i])
       }
       r<-Y - pis
       U<-t(XUs) %*% r
   }

    ##observed statistics
    Ts=rep(NA,length(pow))
    npow = pow
    for (j in 1:length(pow)){
        if (pow[j]<Inf) {
            Ts[j] = sum(U^pow[j])
        } else {
            Ts[j] = max(abs(U))
            npow[j] = 0
        }
    }

    ## residual permutation

    n_pow = length(pow)
    nr_XUs = nrow(XUs)
    nc_XUs = ncol(XUs)
    n_perm = n.perm
    n_r = length(r)
    if(userank)
        r <- jitter(r, amount = 0.0001)

    output <- .C("R_get_pvs",
                 as.double(XUs),
                 as.double(Ts),
                 as.double(npow),
                 as.double(r),
                 as.integer(n_pow),
                 as.integer(nr_XUs),
                 as.integer(nc_XUs),
                 as.integer(n_perm),
                 as.integer(n_r),
                 pvs = as.double( rep(0,n_pow + 1) ),
                 PACKAGE="aSPU")

    pvs <- output$pvs

    Ts <- c(Ts, min(pvs[1:n_pow]) )
    names(Ts) <- c(paste("SPU", pow, sep=""), "aSPU")
    names(pvs) = names(Ts)

    list(Ts = Ts, pvs = pvs)
}
