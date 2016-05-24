## Variance-weighted adptive Sum of powered score (SPUw) test; using permutations to get the p-values (NO adjustment for covariates yet, C version).
##
## It gives the p-values of the SPUw test and aSPUw test based on based on the permutation of residuals.  (NO adjustment for covariates yet, C version)
##
## @param Y phenotype data. It can be disease lables; =0 for controls, =1 for cases.
##     or It can be any quantitative traits. Vector with length n (number of observations)
##
## @param X genotype data; each row for a subject, and each column
##     for an SNP. The value of each element is the # of the copies
##     for an allele.  Matrix with dimension n by k (n : number of observation, k : number of genotype data)
##
## @param cov covariates. Matrix with dimension n by p (n :number of observation, p : number of covariates
##
## @param pow power used in SPU test. Vector of g number of power.
##
## @param n.perm number of permutation
##
## @param userank similar code with the original version if TRUE. by definition if FALSE
##
## @export
## @return Test Statistics and p-values for SPU tests and aSPU test.
##
## @examples
##
## data(exdat)
## out <- aSPUwC(exdat$Y, exdat$X, pow = c(1:8, Inf), n.perm = 1000)
## out
##
## @seealso \code{\link{aSPU}}, \code{\link{aSPUperm2}}, \code{\link{aSPUboot}}, \code{\link{aSPUboot2}}


aSPUwpermC <- function(Y, X, cov = NULL, model=c("gaussian", "binomial"), pow=c(1:8, Inf), n.perm=1000, userank = T){

    model <- match.arg(model)


    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)



    if (is.null(cov)){
        ## NO nuisance parameters:
        XUs <- Xg <- X
        Xbar<-apply(Xg, 2, mean)
        subtract<-function(x, y) { x - y }
        Xgb=t(apply(Xg, 1, subtract, Xbar))

        r=Y-mean(Y)

        U<-as.vector( t(Xg) %*% r)

        ##cov of the score stats:
        if( model == "binomial" ) {
            CovS <- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)
        } else {
            CovS <- var(Y)*(t(Xgb) %*% Xgb)
        }
        
    } else {
        ## with nuisance parameters:
        tdat1<-data.frame(trait=Y, cov)
        fit1<-glm(trait~.,family=model,data=tdat1)
        pis<-fitted.values(fit1)
        XUs<-matrix(0, nrow=n, ncol=k)
        r <- Y - pis
        for(i in 1:k){
            tdat2<-data.frame(X1=X[,i], cov)
            fit2<-glm(X1~.,data=tdat2)
            X1mus<-fitted.values(fit2)
            XUs[, i]<-(X[,i] - X1mus)
        }
        U <- t(XUs) %*% (Y - pis)

        if( model == "binomial" ) {
            CovS <- mean(pis*(1-pis))*(t(XUs) %*% XUs)
        } else {
            CovS <- var(r)*(t(XUs) %*% XUs)
        }        

    }


    Vs<-diag(CovS)
    diagSDs<-ifelse(Vs>1e-20, sqrt(Vs), 1e-10)

   # test stat's:
    Ts<-rep(0, length(pow))
    npow = pow
    for(j in 1:length(pow)){
        if (pow[j] < Inf) {
            Ts[j] = sum((U/diagSDs)^pow[j])
        } else {
            Ts[j] = max(abs(U/diagSDs))
            npow[j] = 0
        }
#     VarTs[j] = var(Upow)
    }

   # permutations:
    n_pow = length(pow)
    nr_XUs = nrow(XUs)
    nc_XUs = ncol(XUs)
    n_perm = n.perm
    n_r = length(r)
    if(userank)
        r <- jitter(r, amount = 0.0001)

    output <- .C("R_get_pvs2",
                 as.double(XUs),
                 as.double(Ts),
                 as.double(npow),
                 as.double(r),
                 as.double(diagSDs),
                 as.integer(n_pow),
                 as.integer(nr_XUs),
                 as.integer(nc_XUs),
                 as.integer(n_perm),
                 as.integer(n_r),
                 pvs = as.double( rep(0,n_pow + 1) ),
                 PACKAGE="aSPU")

    pvs <- output$pvs

    Ts <- c(Ts, min(pvs[1:n_pow]) )
    names(Ts) <- c(paste("SPUw", pow, sep=""), "aSPUw")
    names(pvs) = names(Ts)

    list(Ts = Ts, pvs = pvs)
}
