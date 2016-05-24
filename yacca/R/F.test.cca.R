`F.test.cca` <- function(x, ...)
{
    #Extract basic object information
    if(!("cca"%in%class(x)))
      stop("Object of class cca required.\n")
    s <- length(x$corr)
    p <- length(x$xlab)
    q <- length(x$ylab)
    N <- NROW(x$canvarx)
    k <- 1:s
      
    #Compute statistic and df
    lambda <- sapply(k,function(i){prod(1-x$corr[i:s]^2)})
    r <- (N-s-1)-((abs(p-q)+1)/2)
    Ndf <- (p-k+1)*(q-k+1)
    u <- (Ndf-2)/4
    xx <- ((p-k+1)^2+(q-k+1)^2)-5
    t <- sqrt(((p-k+1)^2*(q-k+1)^2-4)/xx)
    ilambda <- lambda^(1/t)
    Ddf <- (r*t)-(2*u)
    Fstat <- ((1-ilambda)/ilambda)*(Ddf/Ndf)
    pgF <- pf(Fstat,Ndf,Ddf,lower.tail=FALSE)
    
    #Assemble and return the results
    out <- list()
    out$corr <- x$corr
    out$statistic <- Fstat
    out$parameter <- cbind(Ndf,Ddf)
    colnames(out$parameter) <- c("num df", "denom df")
    out$p.value <- pgF
    out$method <- "F test for significance of canonical correlations"
    out$data.name <- names(x$corr)
    class(out)<-c("F.test.cca","htest")
    out
}
