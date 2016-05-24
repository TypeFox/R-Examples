AutoBoot.test <-
function(y,nboot,wild,prob=c(0.025,0.975))
{
    set.seed(12345)
    y <- as.matrix(y)
    test <- Auto.VR(y); LC=test$stat
    
    statmat1 <- matrix(NA, nrow=nboot, ncol=1)
    statmat2 <- matrix(NA, nrow=nboot, ncol=1)
    if (wild == "Normal")
    {
        for (i in 1:nboot)
        {
        ys <- y * rnorm(nrow(y))
        M=Auto.VR(ys)
        statmat1[i,] <- M$stat; statmat2[i,]=M$sum 
        }
    }
    
    if (wild == "Mammen")
    {
        for (i in 1:nboot)
        {
        ys <- y * Mammen(nrow(y))
        M=Auto.VR(ys)
        statmat1[i,] <- M$stat; statmat2[i,]=M$sum 
        }
    }
    
    if (wild == "Rademacher")
    {
        for (i in 1:nboot)
        {
        ys <- y * Rademacher(nrow(y))
        M=Auto.VR(ys)
        statmat1[i,] <- M$stat; statmat2[i,]=M$sum 
        }
    }
     tem <- abs(statmat1) > abs(LC)
     tem[tem == "TRUE"] <- 1
      p <- mean(tem)
      CI1 <- quantile(statmat1,prob);  CI2 <- quantile(statmat2,prob)
return(list(test.stat=LC,VRsum=test$sum,pval=p,CI.stat=CI1,CI.VRsum=CI2))
}
