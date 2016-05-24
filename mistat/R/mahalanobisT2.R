mahalanobisT2 <- function(x, 
                          factor.name, 
                          response.names=names(x)[!names(x) %in% factor.name],
                          conf.level=0.95,
                          compare.to=NA,
                          plot=FALSE){
  x <- x[, c(factor.name, response.names)]
  
  names(x)[1] <- "factor"
  
  D <- as.data.frame(t(aggregate(. ~ factor, data=x, FUN=mean)[, -1]))
  
  Cov <- by(x[, response.names], INDICES=x[,"factor"], cov)
  Cov <- do.call("+", Cov)/length(Cov)
  
  N <- as.data.frame(t(aggregate(. ~ factor, data=x, FUN=length)[, -1]))
  if(any(N != N[1, 1]))
    warning(paste("different number of measures by factor, using n=", N[1, 1], sep=""))
  N <- N[1, 1]
  
  P <- nrow(D)
  
  K <- N^2/(2*N)*(2*N-P-1)/((2*N-2)*P)
  
  Qf <- qf(p=conf.level, df1=P, df2=2*N-P-1)
    
  D2 <- as.matrix(D$V1 - D$V2)
  
  D3 <- rbind(0, t(D2))
  colnames(D3) <- LETTERS[1:ncol(D3)]
  
  Lm <- lm(A ~ -1 + ., data=as.data.frame(D3))
  
  temp <- function(x, lm, d3, d2, cov, k, qf){
    x <- as.data.frame(t(x)) 
    names(x) <- LETTERS[2:max(2, ncol(d3)-1)]
    d4 <- as.numeric(c(predict(lm, newdata=x), x))
    res <- abs(k*((t(d4) - t(d2)) %*% 
                       solve(cov) 
                     %*% (d4 - d2)) - qf)
    return(res)
  }
  
  Par <- optim(par=rep(0, ncol(D3)-1), 
               fn=temp, lm=Lm, d3=D3, d2=D2, cov=Cov, k=K, qf=Qf,
               method=ifelse(ncol(D3)-1 > 1, "Nelder-Mead", "BFGS"), 
               control=list(abstol=0.000001))$par
  
  New <- as.data.frame(t(Par)) 
  names(New) <- LETTERS[2:max(2, ncol(D3)-1)]
  
  New <- as.numeric(c(predict(Lm, newdata=New), New))
  res <- cbind(New, D2, D2-New+D2)
  
  res <- list(coord = t(res[,order(res[1,])]))
  
  colnames(res$coord) <- response.names
  rownames(res$coord) <- c("LCR", "Center", "UCR")
  
  res$mahalanobis <- apply(res$coord, 1, FUN=function(x, cov){sqrt(t(x) %*% solve(cov) %*% x)}, cov=Cov)
  names(res$mahalanobis) <- c("LCR", "Center", "UCR")
  
  if(!any(is.na(compare.to))){
    res$mahalanobis.compare <- sqrt(t(compare.to) %*% solve(Cov) %*% compare.to)
  } else {
    res$mahalanobis.compare <- as.numeric(NA)
  }
  
  if(plot && ncol(res$coord) <= 2){
    OldPar <- par("mar")
    layout(matrix(1:2), heights=c(2,1))
    par(mar=c(4, 4, 1, 2.1))
    temp <- rbind(0, res$coord)
    plot(temp, main="")
    lines(temp[c(2, nrow(temp)), ], col=2)
    lines(temp[c(1, 2), ])
    points(res$coord, pch=20)
    grid()
    par(mar=c(1, 4, 4, 2.1))
    plot(c(0:2, 0), 
         c(res$mahalanobis, res$mahalanobis.compare), 
         type="n", axes=FALSE, ylab=expression("Mahalanobi T^2"), xlab="")
    axis(side=2)
    rect(0.5, min(res$mahalanobis), 1.5, max(res$mahalanobis), col=2)
    lines(x=c(0.5, 1.5), rep(res$mahalanobis[2], 2), col=1, lwd=3)
    if(!any(is.na(compare.to))){
      lines(x=c(0.25, 1.6), rep(res$mahalanobis.compare, 2), col=3, lwd=3)
      text(x=1.75, y=res$mahalanobis.compare, labels="Comparison")
    }
    par(OldPar)
    layout(1)
  }
  return(res)  
}