## Sum of powered score (SPU) test (boot, version 2, matrix used in permutation)
##
## It gives the p-values of the SPU test and aSPU test based on the parametric bootstrap. (This is version 2, matrix version is faster but if it doesn't work, we should use version 1, vector version)
##
## @param Y phenotype data. It can be disease lables; =0 for controls, =1 for cases.
## or it can be any quantitative traits. Vector with length n (number of observations)
##
## @param X genotype data; each row for a subject, and each column
##     for an SNP. The value of each element is the # of the copies
##     for an allele. Matrix with dimension n by k (n : number of observation, k : number of genotype data)
##
## @param cov covariates. Matrix with dimension n by p (n :number of observation, p : number of covariates)
##
## @param model Use "gaussian" for quantitative trait, and use "binomial" for binary trait.
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
## out <- aSPUboot2(exdat$Y, exdat$X, cov = NULL,
##                model = "binomial", pow = c(1:8, Inf), n.perm = 1000)
## out
##
## @seealso \code{\link{aSPU}}, \code{\link{aSPUperm}}, \code{\link{aSPUperm2}}, \code{\link{aSPUboot}}



aSPUboot2 <- function(Y, X, cov=NULL, model=c("gaussian", "binomial"), pow=c(1:8, Inf), n.perm=1000){

    model <- match.arg(model)

    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)

    if (is.null(cov)){
        ## NO nuisance parameters:
        Xg <- X
        U <- t(Xg) %*% (Y-mean(Y))
        yresids <- Y-mean(Y)
#        sigma0 = sqrt(sum(yresids^2)/(n-1))
        yfits <- rep(mean(Y), n)

    } else {
        ## with nuisance parameters:
        tdat1 <- data.frame(trait=Y, cov)
        fit1 <- glm(trait~., family = model, data=tdat1)
        yfits <- fitted.values(fit1)
        yresids <- Y - yfits
 #       fit1res1<-summary(fit1)
 #       sigma0<-sqrt(fit1res1$dispersion)

        Us <- XUs <- matrix(0, nrow=n, ncol=k)
        Xmus = X
        for(i in 1:k){
            tdat2 <- data.frame(X1=X[,i], cov)
            fit2 <- glm(X1~., data=tdat2)
            Xmus[,i] <- fitted.values(fit2)
            XUs[, i] <- (X[,i] - Xmus[,i])
        }
        U <- t(XUs) %*% (Y - yfits)
    }
                                        # test stat's:
    Ts <- rep(0, length(pow))
    for(j in 1:length(pow)){
        if (pow[j] < Inf)
            Ts[j] = sum(U^pow[j]) else Ts[j] = max(abs(U))
    }

                                        # bootstrap:
    T0s = matrix(0, nrow=n.perm, ncol=length(pow))
    Y0 = Y
    for(b in 1:n.perm){
        if (is.null(cov)) {
            Y0 <- sample(Y, length(Y))
#########Null score vector:
            U0 <- t(Xg) %*% (Y0-mean(Y0))
        }
        else{
            ## with nuisance parameters:
            if ( model == "gaussian") {

                Y0 <- yfits + sample(yresids, n, replace = F )
                tdat0<-data.frame(trait=Y0, cov)
                fit0<-glm(trait~., data=tdat0)
                yfits0<-fitted.values(fit0)
                U0<-t(XUs) %*% (Y0 - yfits0)
            } else {
                ## with nuisance parameters:
                for(i in 1:n) Y0[i] <- sample(c(1,0), 1, prob=c(yfits[i], 1-yfits[i]) )
                tdat0<-data.frame(trait=Y0, cov)
                fit0<-glm(trait~., family=model, data=tdat0)
                yfits0<-fitted.values(fit0)
                U0<-t(XUs) %*% (Y0 - yfits0)
            }

        }

                                        # test stat's:
        for(j in 1:length(pow))
            if (pow[j] < Inf)
                T0s[b, j] = sum(U0^pow[j]) else T0s[b, j] = max(abs(U0))

    }

                                        # bootstrap-based p-values:
                                        #pPerm0 <- apply( matrix( rep(abs(Ts),n.perm), nrow = n.perm, byrow = T) < abs(T0s), 2, mean)

    pPerm0 = rep(NA,length(pow))
    for ( j in 1:length(pow))
    {
        pPerm0[j] = round( sum(abs(Ts[j])<=abs(T0s[,j])) / n.perm, digits = 8)
        P0s = ( ( n.perm - rank( abs(T0s[,j]) ) ) + 1 ) / (n.perm)
        if (j == 1 ) minp0  = P0s else minp0[which(minp0>P0s)] = P0s[which(minp0>P0s)]
    }

    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    pvs <- c(pPerm0, Paspu)

    Ts <- c(Ts, min(pPerm0))
    names(Ts) <- c(paste("SPU", pow, sep=""), "aSPU")
    names(pvs) = names(Ts)

    list(Ts = Ts, pvs = pvs)

}
