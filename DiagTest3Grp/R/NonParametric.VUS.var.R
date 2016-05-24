NonParametric.VUS.var <-
function(x,y,z,alpha=0.05,NBOOT=100,FisherZ=FALSE)
  {
    ######################################################################################################################################################
    ################This function provide variance on nonparametric estimate of VUS by bootstrapping sampling##############
    ########
    ####1. Input:
    ###(1)x,y,z:vectors are the test marker measurements for the 3 ordinal group (D-,D0,D+)
    ###(2)alpha: provide  bootstrap basic (1-alpha)% confidence interval, alpha=0.05 gives (5%,97.5%) CI
    ###(3)NBOOT: number of bootstrapping
    ###(4)FisherZ: default=FALSE, if TRUE, calculate the variance and quantile of VUS on the logit scale
    ####2. Output: (1) variance on nonparametric VUS estimate from bootstrapping (2) bootstrap basic (1-alpha)% quantile interval
    ######################################################################################################################################################
    
    #cat("It may take a minute to obtain bootstrapping variance\n")  
    boot.VUS <- sapply(1:NBOOT,function(jj)
                       {

                         new.x <- bootSample(x,jj)##use jj as seed for random number generation
                         new.y <- bootSample(y,jj+1)##use jj+1 as seed for random number generation
                         new.z <- bootSample(z,jj+2)##use jj+2 as seed for random number generation
                         res0 <- NonParametric.VUS(x=new.x,y=new.y,z=new.z,FisherZ=FisherZ)
                         res0$estimate
                       })
    ###caculate variance
    var0 <- var(boot.VUS,na.rm=T)
    prob0 <- alpha/2
    ###CI
    lower <- quantile(boot.VUS,prob=prob0)
    upper <- quantile(boot.VUS,prob=1-prob0)

    CI <- c(lower,upper)
    names(CI) <- c(paste(alpha/2*100,"%",sep=""),paste(100-alpha/2*100,"%",sep=""))
    
    return(list(variance=var0,CI=CI))
  }

