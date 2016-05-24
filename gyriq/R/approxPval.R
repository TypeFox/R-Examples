approxPval <- function(score, K, asv=NULL, method="davies", starResid=NULL)
## Approximates a p-value when the test statistic is a quadratic form following
## a mixture of chi-squared variables.
## 
## 'score': vector composing the quadratic form
## 'K': matrix composing the quadratic form
## 'asv': (default=NULL) number of approximate eigenvalues to be estimated for 
##      the matrix if the implicitly-restarted Lanczos bidiagonalization is used
## 'method': (default="davies") procedure used to obtain the p-value
## 'starResid': (default=NULL) matrix of permuted residuals used to obtain the
##      p-value if a permutation procedure is employed
{
    if (method == "rspOrd") {
        return(list(aPvalResampOrd(q=score, K=K, starResid=starResid), 
                    KeigVal=NULL))
    }
    
    if (method == "rspMom") {
        return(list(aPvalResampMom(q=score, K=K, starResid=starResid), 
                    KeigVal=NULL))
    }

    if (method == "davies") {
        if (is.null(asv) == TRUE) {
            KeigVal <- eigen(K)$values
        } else {
            KeigVal <- irlba::irlba(K, nu=asv, nv=asv)$d
        }

        eigIndex <- seq(along.with=KeigVal)[Re(KeigVal) > 0 & Im(KeigVal) == 0]
        muVal <- as.numeric(KeigVal[eigIndex])

        return(list(CompQuadForm::davies(q=score, lambda=muVal)$Qq, KeigVal))
    }
}