sd_trim <- function(x,trim=0.2, const=TRUE){
      # trimmed sd, where x is a matrix (column-wise)
      x <- as.matrix(x)
      if (const){
        if (trim==0.1){const <- 0.7892}
        else if (trim==0.2){const <- 0.6615}
        else {warning("Did you specify the correct consistency constant for trimming?")}
      }
      else{const <- 1}
      m <- apply(x,2,mean,trim)
      res <- x-rep(1,nrow(x))%*%t(m)
      qu <- apply(abs(res),2,quantile,1-trim)
      sdtrim <- apply(matrix(res[t(abs(t(res))<=qu)]^2,ncol=ncol(x),byrow=FALSE),2,sum)
      sdtrim <- sqrt(sdtrim/((nrow(x)*(1-trim)-1)))/const
      return(sdtrim)
}

