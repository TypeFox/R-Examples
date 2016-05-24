#
mauchly.test <- function(X,location, n, p)
    {
    if (location=="est") S <- cov(X) else S <- covOrigin(X)
    L <- ((det(S) / (sum(diag(S))/p)^p ))^(n/2)
    names(L) <- "L"
    df <- (p+2)*(p-1)/2
    names(df)<-"df"
    p.value <- 1-pchisq(-2*log(L), df)
    METHOD<-"Mauchly test for sphericity"
    list(statistic=L, p.value = p.value, parameter = df, method = METHOD)
    }

SSCov.test <- function(X, n, p)
    {
    cpp <- SpatialNP:::Cpp(p)
    METHOD = "Test of sphericity based on TCOV"
    TCOV <- SSCov(X)
    vTCOV <- as.vector(TCOV)
    cross <- SpatialNP:::Q2internal(X)
    cov.cTCOV <- 4 * (matrix(cross[-(1:p^2)], ncol = p^2) -  tcrossprod(vTCOV)) /n
    test.statistic <- as.vector(t(cpp %*% vTCOV) %*% SpatialNP:::gen.inv(cov.cTCOV) %*% (cpp %*% vTCOV))
    names(test.statistic) <- "Q2"
    df <- (p+2)*(p-1)/2
    names(df)<-"df"
    p.value <- 1-pchisq(test.statistic, df)
    METHOD<-"Test for sphericity based on TCOV"
    list(statistic=test.statistic, p.value = p.value, parameter = df, method = METHOD)
    }

SCov.test <- function(X,location, n, p)
    {
    if (location=="est") T <- spatial.sign(X,center=TRUE,shape=FALSE) else T <- spatial.sign(X,center=FALSE,shape=FALSE)
    vSCOV <- as.vector(crossprod(T)/n)
    cpp<- SpatialNP:::Cpp(p)
    Tau <- (2 / (p*(p+2))) 
    Q2 <- sum((cpp %*% vSCOV)^2)
    test.statistic <- n*Q2/Tau
    df <- (p+2)*(p-1)/2
    p.value <- 1-pchisq(test.statistic, df)
    names(df)<-"df"
    names(test.statistic) <- "Q2"
    METHOD<-"Test for sphericity based on UCOV"
    list(statistic=test.statistic, p.value = p.value, parameter = df, method = METHOD)
    }
    
    
