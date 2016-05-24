## Adaptive Sum of powered score (SPU) tests (SPU and aSPU) (simulation, version 1, vector used in permutation)
##
## It gives the p-values of the SPU tests and aSPU test based on the simulations of U from the null distribution. (This is version 1, matrix version is faster but if it doesn't work, we should use version 1, vector version)
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
## @export
## @return Test Statistics and p-values for SPU tests and aSPU test.
##
## @examples
##
## data(exdat)
## out <- aSPUsim1(exdat$Y, exdat$X, cov = NULL, model = "binomial",
##                   pow = c(1:8, Inf), n.perm = 1000)
## out
##
## @seealso \code{\link{aSPU}}, \code{\link{aSPUperm2}}, \code{\link{aSPUboot}}, \code{\link{aSPUboot2}}


aSPUsim1 <- function(Y, X, cov = NULL, model=c("gaussian","binomial"), pow=c(1:8, Inf), n.perm=1000){

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

        Xbar<-apply(Xg, 2, mean)
        Xgb<-Xg
        for(i in 1:nrow(Xg))
            Xgb[i,]<-Xg[i,]-Xbar

	if( model == "binomial" ) {
            CovS <- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)
	} else {
            CovS <- var(Y)*(t(Xgb) %*% Xgb)
	}
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
        
        if( model == "binomial" ) {
            CovS <- mean(pis*(1-pis))*(t(Xgb) %*% Xgb)
        } else {
            CovS <- var(r)*(t(Xgb) %*% Xgb)
        }
        
#        CovS<-matrix(0, nrow=k, ncol=k)
#        for(i in 1:n)
#            CovS<-CovS + XUs[i,] %*% t(XUs[i,])
    }

    svd.CovS<-svd(CovS)
    CovSsqrt<-svd.CovS$u %*% diag(sqrt(svd.CovS$d))


    ##observed statistics
    Ts=rep(NA,length(pow))
    for (j in 1:length(pow)){
        if (pow[j]<Inf) Ts[j] = sum(U^pow[j]) else Ts[j] = max(abs(U))
    }
    ## cat("statistic calculated","\n")


    ## simulation based
    pPerm0 = rep(NA,length(pow))
    T0s = numeric(n.perm)
    s <- sample(1:10^5,1)

    for (j in 1:length(pow)){
        set.seed(s) # to ensure the same samples are drawn for each pow
        for (b in 1:n.perm){
            U00<-rnorm(k, 0, 1)
            U0<-CovSsqrt %*% U00
#            r0 <- sample(r, length(r))
#            U0 <- as.vector(t(XUs) %*% r0)
            if (pow[j] < Inf){ T0s[b] = round(sum( U0^pow[j]), digits = 8) }
            if (pow[j] == Inf) {T0s[b] = round( max(abs(U0)), digits = 8) }
	}

        pPerm0[j] = round( sum(abs(Ts[j])<=abs(T0s)) / n.perm, digits = 8)
        P0s = ( (n.perm-rank(abs(T0s))) + 1 )/(n.perm)
        if (j==1) minp0=P0s else minp0[which(minp0>P0s)]=P0s[which(minp0>P0s)]
    }

#    cat("P0s caculated","\n")
    Paspu<-(sum(minp0<=min(pPerm0))+1)/(n.perm+1)
    pvs <- c(pPerm0, Paspu)

    Ts <- c(Ts, min(pPerm0))
    names(Ts) <- c(paste("SPU", pow, sep=""), "aSPU")
    names(pvs) = names(Ts)

    list(Ts = Ts, pvs = pvs)
}

