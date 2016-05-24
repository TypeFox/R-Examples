Boot.test <-
function(y,kvec,nboot,wild,prob=c(0.025,0.975))
{
    set.seed(12345)
    y <- as.matrix(y)
    LC <- FastLMCD(y,kvec)
    
    statmat <- matrix(NA, nrow=nboot, ncol=length(kvec)+1)
    
    if (wild == "Normal")
    {
        stat <- matrix(c(LC$M2,LC$CD2))
        for (i in 1:nboot)
        {
        ys <- y * rnorm(nrow(y))
        LCs <- FastLMCD(ys,kvec)
        statmat[i,] <- c(LCs$M2,LCs$CD2)
        }
    }

    if (wild == "No")
    {
        stat <- matrix(c(LC$M1,LC$CD1))
        for (i in 1:nboot)
        {
        index <- as.integer(runif(nrow(y), min=1, max=nrow(y)))
        ys <- as.matrix(y[index])
        LCs <- FastLMCD(ys,kvec)
        statmat[i,] <- c(LCs$M1,LCs$CD1)
        }
    }
    
    if (wild == "Mammen")
    {
        stat <- matrix(c(LC$M2,LC$CD2))
        for (i in 1:nboot)
        {
        ys <- y * Mammen(nrow(y))
        LCs <- FastLMCD(ys,kvec)
        statmat[i,] <- c(LCs$M2,LCs$CD2)
        } 
    
    }
    
    if (wild == "Rademacher")
    {
        stat <- matrix(c(LC$M2,LC$CD2))
        for (i in 1:nboot)
        {
        ys <- y * Rademacher(nrow(y))
        LCs <- FastLMCD(ys,kvec)
        statmat[i,] <- c(LCs$M2,LCs$CD2)
        } 
    
    }
    
    p <- matrix(NA,nrow = ncol(statmat), ncol=1)
    CI <- matrix(NA,nrow = ncol(statmat), ncol=length(prob))
    for (i in 1:ncol(statmat))
        {
        tem <- abs(statmat[,i]) > abs(stat[i])
        tem[tem == "TRUE"] <- 1
        p[i] <- mean(tem) 
        CI[i,] <- quantile(statmat[,i],prob)
        }
   CI <- CI[1:length(kvec),]
   colnames(CI) <- paste(100*prob,"%",sep="")
   rownames(CI) <- paste("k=",kvec,sep="")
         
return(list(Holding.Period=kvec,LM.pval=as.numeric(p[1:length(kvec)]),CD.pval=as.numeric(p[length(kvec)+1]),CI=CI))
}
