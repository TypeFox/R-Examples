Subsample.test <-
function(y,kvec)
{
    y <- as.matrix(y)
    n <- nrow(y)
    b1 <- as.integer(2.5*n^(0.3))
    b2 <- as.integer(3.5*n^(0.6))
    term <- as.integer( (b2-b1)/7 ) 
    b1vec <- as.matrix(seq(b1,b2,term)[2:7])

    p <- matrix(NA,nrow = nrow(b1vec), ncol=1)

    for (i in 1:nrow(b1vec))    
    {   
        mv <- WK_stat2(y,kvec)
        b1 <- b1vec[i]
        mvsamp <- matrix(NA,nrow=(n-b1+1),ncol=1) 
        index <- 1:b1
        for (j in 1:(n-b1+1))
            {
            xsub <- as.matrix(y[index])
            mvsamp[j] <- WK_stat2(xsub,kvec)
            index <- index+1
            }
   
   tem <- mvsamp > mv
   tem[tem == "TRUE"] <- 1
   p[i] <- mean(tem) 
   }   
   rownames(p) <- paste("bl=",b1vec,sep="")
   colnames(p) <- c("pval")
   return(list(Holding.Period=kvec,Pval=p,Block.length=as.numeric(b1vec)))
}
