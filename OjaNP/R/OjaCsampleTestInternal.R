oja2sampleSignTest <- function(X, Y, p, method, n.simu, center, ...)
    {
    n1 <- nrow(X)
    n2 <- nrow(Y)
    N <- n1+n2
    k <- ncol(X)
    
    SIGNS <- ojaSign(rbind(X,Y),p=p, center=center, ...)
    a.i <- rep(c(1,-1), c(n1,n2))
    T.N <- colSums(a.i*SIGNS)
    V.N <- 4*n1*n2 / (N*(N-1)) * crossprod(SIGNS)
    V.N.inv <- solve(V.N)
    statistic <- as.numeric(t(T.N) %*% V.N.inv %*% T.N) 
    names(statistic) <- "Q.S"
    
    simu.f <- function(index,a.i,SIGNS, V.inv)
        {
        signs <- a.i[index]
        T.N.j <- colSums(signs*SIGNS)
        Q.S.j <- t(T.N.j) %*% V.inv %*% T.N.j
        as.numeric(Q.S.j) 
        }
    
    
    if (method=="approximation") {
        parameter<-k
        names(parameter)<-"df"   
        p.value<-1-pchisq(statistic,k)
        } else{
        stats <- replicate(n.simu, simu.f(index = sample(1:N,N), a.i=a.i, SIGNS=SIGNS, V.inv=V.N.inv))
        parameter <- n.simu
        names(parameter) <- "permutations"
        p.value <- mean(statistic < stats)
        }
    
    METHOD <- "OJA C SAMPLE SIGN TEST"
    
    list(statistic=statistic, p.value=p.value, parameter=parameter, method=METHOD)
    
    }


oja2sampleRankTest <- function(X, Y, p, method, n.simu, ...)
    {
    n1 <- nrow(X)
    n2 <- nrow(Y)
    N <- n1+n2
    k <- ncol(X)
    
    RANKS <- ojaRank(rbind(X,Y), p=p, ...)
    LAMBDA <- n2/N
    a.i <- rep(c(-LAMBDA, 1-LAMBDA), c(n1,n2))
    T.N <- colSums(a.i*RANKS)
    B.N <-  crossprod(RANKS) / (N-1)
    B.N.inv <- solve(B.N)
    statistic.1 <- as.numeric(t(T.N) %*% B.N.inv %*% T.N) 
    statistic <- statistic.1 / (N * LAMBDA * (1-LAMBDA))
    names(statistic) <- "Q.R"
    
    simu.f <- function(index,a.i,RANKS, B.inv)
        {
        signs <- a.i[index]
        T.N.j <- colSums(signs*RANKS)
        Q.R.j <- t(T.N.j) %*% B.inv %*% T.N.j
        as.numeric(Q.R.j) 
        }
    
    
    if (method=="approximation") {
        parameter<-k
        names(parameter)<-"df"   
        p.value<-1-pchisq(statistic,k)
        } else{
        stats <- replicate(n.simu, simu.f(index = sample(1:N,N), a.i=a.i, RANKS=RANKS, B.inv=B.N.inv))
        parameter <- n.simu
        names(parameter) <- "permutations"
        p.value <- mean(statistic.1 < stats)
        }
    
    METHOD <- "OJA C SAMPLE RANK TEST"
    
    list(statistic=statistic, p.value=p.value, parameter=parameter, method=METHOD)
    
    }
