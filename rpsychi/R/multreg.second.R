multreg.second <- function(formula, corr, n, m = NULL, sd=NULL, sig.level=0.05, digits=3){

    ##setting
    depname <- unlist(strsplit(as.character(formula), " ")[2])
    x <- unlist(strsplit(as.character(formula), " ")[3])
    indname <- x[!is.element(x, "+")]
    x <- c(depname, indname)
    corr <- corr[x, x]
    sd   <- sd[x]
    m    <- m[x]


    ##(a) correlation
    corr.lower <- tanh(atanh(corr) + qnorm(sig.level/2, lower.tail=TRUE) / sqrt(n - 3))
    corr.upper <- tanh(atanh(corr) + qnorm(sig.level/2, lower.tail=FALSE) / sqrt(n - 3))
    corr.conf  <- corr.lower
    corr.conf[upper.tri(corr.conf)] <- corr.upper[upper.tri(corr.upper)]


    ##(b) partial correlation cor
    p <- ncol(corr)
    K <- solve(corr)
    a <- 1/sqrt(diag(K))
    K <- K * outer(a, a)
    partial.corr <- 2 * diag(p) - K
    dimnames(partial.corr) <- dimnames(partial.corr)

    cor.mat <- corr
    cor.mat[upper.tri(cor.mat)] <- partial.corr[upper.tri(partial.corr)]
    
      
    
    ##(c) standardized coefficients  
    num <- which(depname==colnames(corr)) #check dependent and independent variables
    rxy  <- corr[,num][-num]              #correlations between dependent and independent variables
    Rind <- corr[-num,-num]               #correlations among independent variables
    bs <- solve(Rind) %*% rxy             #standardized coefficients  
    R.sq  <- as.vector(t(as.matrix(rxy)) %*% solve(Rind) %*% as.matrix(rxy))
    
    k <- nrow(bs)
    bs.sem <- numeric(k)
    
    for(i in 1:k){
      xname <- rownames(bs)[i]
      xdel <- setdiff(rownames(bs), xname)
            
      num <- which(xname==colnames(Rind))   #check dependent and independent variables
      rxy  <- Rind[,num][-num]              #correlations between dependent and independent variables
      Rind2 <- Rind[-num,-num]              #correlations among independent variables
      Ri <- as.vector(t(as.matrix(rxy)) %*% solve(Rind2) %*% as.matrix(rxy))
      bs.sem[i] <- sqrt((1 - R.sq)/(n- k -1)) * sqrt(1/(1 - Ri))           
    }

    standardized.estimates <- data.frame(matrix(NA, nrow=k, ncol=4))
    rownames(standardized.estimates) <- rownames(bs)
    colnames(standardized.estimates) <- c("estimates", "lower", "upper", "std")
    standardized.estimates[,1] <- as.numeric(bs)
    standardized.estimates[,4] <- bs.sem
    standardized.estimates[,2] <- standardized.estimates$estimates + qnorm(sig.level/2) * standardized.estimates$std
    standardized.estimates[,3] <- standardized.estimates$estimates + qnorm(sig.level/2, lower.tail=FALSE) * standardized.estimates$std    
    
    
    ##(d) raw coefficients
    if(!is.null(sd)){
      b     <- sd[depname]/sd[rownames(bs)] * standardized.estimates$estimates
      b.sem <- sd[depname]/sd[rownames(bs)] * standardized.estimates$std
      b.lower <- b + qt(sig.level/2, n - k - 1) * b.sem
      b.upper <- b + qt(sig.level/2, n - k - 1, lower.tail=FALSE) * b.sem
      
      Intercept <- m[depname] - sum(b * m[names(b)])
      raw.estimates <- data.frame(matrix(NA, nrow=k+1, ncol=4))
      rownames(raw.estimates) <- c("Intercept", names(b))
      colnames(raw.estimates) <- c("estimates", "lower", "upper", "std")      
      raw.estimates[,1] <- c(Intercept, b)
      raw.estimates[,4] <- c(NA, bs.sem)
      raw.estimates[,2] <- c(NA, b.lower)
      raw.estimates[,3] <- c(NA, b.upper)
    }
        
    
    ##(e) omnibus effects    
    u <- length(indname)
    nu <- n - u - 1
    f.sq <- R.sq/(1-R.sq)
    f.value <- f.sq * (nu / u)

    delta.lower <- try(FNONCT(f.value, u, nu, prob=1-sig.level/2))
    delta.upper <- FNONCT(f.value, u, nu, prob=sig.level/2)
    if(is.character(delta.lower)){
      delta.lower <- 0
    }
    R.sq.lower  <- delta.lower / (delta.lower+u+nu+1)
    R.sq.upper  <- delta.upper / (delta.upper+u+nu+1)
    omnibus.es <- c(Rsq=R.sq, lower=R.sq.lower, upper= R.sq.upper)

    
    
    ##(f) power
    criterion.power <- c(
            small=power.multi(n=n, n.ind=u, delta=0.02, sig.level=sig.level), 
            medium= power.multi(n=n, n.ind=u, delta=0.15, sig.level=sig.level), 
            large = power.multi(n=n, n.ind=u, delta=0.35, sig.level=sig.level)
            )
    
       
  
    ##output
    output <- list(corr.partial.corr=cor.mat, corr.confidence = corr.conf, omnibus.es=omnibus.es, standardized.estimates=standardized.estimates, power=criterion.power)
    if(!is.null(sd)){
    output <- list(corr.partial.corr=cor.mat, corr.confidence = corr.conf, omnibus.es=omnibus.es, raw.estimates=raw.estimates,standardized.estimates=standardized.estimates, power=criterion.power)    
    }
    
    output <- sapply(output, round, digits)
    return(output)        
}