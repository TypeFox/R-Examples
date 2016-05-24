
print.DiscML <- function (x, digits = 6, ...) 
{
  if(identical(x$individualrates,TRUE) && x$NROW >1)
  { 
    #  y <- (x$total) 
    #   y[!is.na( suppressWarnings(matrix( as.numeric(( as.matrix(y)) ), nrow(y),ncol(y))))] <-formatC(round(as.numeric(y[!is.na(suppressWarnings(matrix(as.numeric(( as.matrix(y)) ), nrow(y),ncol(y))))]),4), format= 'f', digits=4)
    #   y[,1] <- round(1:length(y[,1]),0)
    options(scipen= 5)
    cat("\n")
    cat("\n  DiscML (Maximum Likelihood Method for Discrete Character States)\n\n")
    cat("Call: ")
    print(x$call)
    cat("\n")
    ratemat <- x$index.matrix
    dimnames(ratemat)[1:2] <- dimnames(x$lik.anc)[2]
    cat("Rate index matrix:\n")
    colnames(ratemat)<-x$lvls
    rownames(ratemat)<-x$lvls
    print(ratemat, na.print = ".", quote=FALSE)
    cat("\nResult summary (type '...$total' to get this data): \n\n")
    qq<- matrix(0, x$NROW , length(colnames(x$dlist[[1]])))
    for(i in 1:x$NROW)
      qq[i, ] <- as.matrix(x$dlist[[i]])
    qq<- rbind(colnames(x$dlist[[1]]), qq)
    for(i in 1:x$NROW)
    {
      tt<- x$dlist[[i]]
      dd <-(as.matrix(tt))
      colnames(dd) <- NULL    
      dd[ nchar(dd) > 10] <- formatC( as.numeric(dd[ nchar(dd) > 10]) , format = 'e', digits= digits )   
      qq[i+1, ]<-dd
    }
    for(i in 1:(x$NROW+1))
    {
      dd<- qq[i, ]
    for(j in 1:length(dd))
    {
      dd[j] <- paste(dd[j] , paste( rep(" ", max( nchar(qq[ , j]) ) - nchar(dd[j]) ) ,collapse=""   ) ,   collapse="")
    }
    
    cat((dd),"\n")
    }
    cat("\nTotal time elapsed in running this function (type '...$time' to get this data): \n")
    cat(  paste( round(x$time, 2), "seconds"),"\n")
    if(FALSE)
    {
      if (!is.null(x$loglik)) 
      {if(x$NROW >1 )
      {
        cat("Log Likelighood (type '...$loglik' to get this data):\n\n")
        print(round(x$loglik, 8))
        cat("\n")
      }
      else
      {
        cat("Log Likelighood (type '...$loglik' to get this data): ")
        print(round(x$loglik, 8))
        cat("\n\n")
      }
      }
      ratemat <- x$index.matrix
      dimnames(ratemat)[1:2] <- dimnames(x$lik.anc)[2]
      cat("Rate index matrix:\n")
      colnames(ratemat)<-x$lvls
      rownames(ratemat)<-x$lvls
      print(ratemat, na.print = ".", quote=FALSE)
      cat("\n")
      npar <- length(x$rates[,1])
      
      if(identical(x$ismu,TRUE))
        cat("Standardized matrix values:\n ")  
      else
        cat("Rate parameter estimates (type '...$rates' to get this table):\n ")
      print(  round(x$rates, digits), row.names= FALSE)
      
      if(identical(x$ismu, TRUE))
      {
        cat("\nThe values of mu\'s estimated (Type ...$mu to get this table):\n")           
        x$mu <- cbind(x$mu[,1, drop =FALSE] , round(x$mu[,-1, drop =FALSE], digits ))
        print(x$mu, row.names = FALSE)
      }
      if(identical(x$isrprob, TRUE))
        cat("\nPrior root probability estimated (type '...$rprob' to get this table): \n")
      else
        cat("\nPrior root probability used (type '...$rprob' to get this table): \n")     
      rprob <- x$rprob
      if(is.data.frame(rprob))
      {
        temp <- sapply(rprob ,is.numeric)
        rprob[ ,temp] <- round(rprob[ ,temp] ,digits)
      }
      print(rprob, row.names = FALSE, quote =FALSE)
      
      #if(identical(x$isrprob, TRUE))
      if(FALSE)
      {
        cat("\nParametric variables estimates (type '...$srprob' to get this table): \n")
        srprob <- x$srprob
        if(is.data.frame(srprob))
        {
          temp <- sapply(srprob ,is.numeric)
          srprob[ ,temp] <- round(srprob[ ,temp] ,digits)
        }
        print(srprob, row.names = FALSE, quote =FALSE)
      }
      cat("\n\n")
      if(!identical(x$isalpha,FALSE)){
        cat("The value of alpha parameter of gamma estimated (or used (type '...$alpha' to get this table):\n")
        print(round(x$alpha,digits))
      }
    }
  }
  if(identical(x$individualrates, FALSE)){ 
    cat("\n  DiscML (Maximum Likelihood Method for Discrete Character States)\n\n")
    cat("Call: ")
    print(x$call)
    cat("\n")
    if (!is.null(x$loglik)) 
      cat("Log Likelighood (type '...$loglik' to get this data): ",x$loglik,"\n")
    cat("\n")
    if (!is.null(x$resloglik)) 
      cat("    Residual log-likelihood:", x$resloglik, "\n\n")
    ratemat <- x$index.matrix
    if (is.null(ratemat)) {
      class(x) <- NULL
      x$resloglik <- x$loglik <- x$call <- NULL
      print(x)
    }
    else {
      dimnames(ratemat)[1:2] <- dimnames(x$lik.anc)[2]
      cat("Rate index matrix:\n")
      colnames(ratemat)<-x$lvls
      rownames(ratemat)<-x$lvls
      print(ratemat, na.print = ".", quote=FALSE)
      cat("\n")    
      if(identical(x$ismu,TRUE))
        cat("Standardized matrix values:\n ")  
      else
        cat("Rate parameter estimates (type '...$rates' to get this table):\n ")
      print( round( x$rates,digits), row.names= FALSE)
      
      if(identical(x$ismu, TRUE))
      {
        cat("\nThe values of mu\'s estimated (Type ...$mu to get this table):\n")           
        x$mu <- cbind(x$mu[,1, drop =FALSE] , round(x$mu[,-1, drop =FALSE], digits ))
        print(x$mu, row.names = FALSE)
      } 
      
      if(identical(x$isrprob, TRUE))
        cat("\nPrior root probability estimated (type '...$rprob' to get this table): \n")
      else
        cat("\nPrior root probability used (type '...$rprob' to get this table): \n")     
      rprob <- x$rprob
      if(is.data.frame(rprob))
      {
        temp <- sapply(rprob ,is.numeric)
        rprob[ ,temp] <- round(rprob[ ,temp] ,digits)
      }
      print( rprob , row.names = FALSE)
      
      #if(identical(x$isrprob, TRUE))
      if(FALSE)
      {
        cat("\nParametric variables estimates (type '...$srprob' to get this table): \n")
        srprob <- x$srprob
        if(is.data.frame(srprob))
        {
          temp <- sapply(srprob ,is.numeric)
          srprob[ ,temp] <- round(srprob[ ,temp] ,digits)
        }
        print(srprob, row.names = FALSE, quote =FALSE)
      }
    }
    if(!identical(x$isalpha,FALSE)){
      if(identical(x$isalpha,TRUE)) 
        cat("\nalpha parameter estimate for discrete-gamma estimate (type '...$alpha' to get this table):\n")
      else
        cat("\nalpha parameter used for discrete-gamma estimate (type '...$alpha' to get this table):\n")
      print(round(x$alpha,digits))
    }
    if (!is.null(x$lik.anc)) {
      writeLines("\nScaled likelihoods at the root for each site 
          (type '...$lik.anc' to get this table):")
      print(round(Re(x$lik.anc), 7) )
    }
    cat("\nTotal time elapsed in running this function (type '...$time' to get this data): \n")
    cat(paste(round(x$time,2), "seconds"), "\n") 
  }
}
