`poly3est` <-
function(time, status, tmax, method="BP", k=NULL)
{

if(is.null(k))
 {k<-3}

# default weight is one

weight <- rep(1,length(status))

    n<-length(status)
    Y<-sum(status)

# replace all weights where no tumor is present with the poly3weight (t.m/tmax)^3

for(i in which(!status))
 {weight[i] <- (time[i]/tmax)^k}

    nadj<-sum(weight)
    estimate<-Y/nadj 

if (method=="BP")
    {
    ##Varianzberechnung Bailer Portier
    estp<-estimate
    varp<-estp*(1-estp)/nadj 
    varcor<-varp

    if(Y==0)
    {
    estcor<-0.5/nadj
    varcor<-estcor*(1-estcor)/nadj
    }

    if(Y==n)
    {
    estcor<- (n-0.5)/nadj
    varcor<-estcor*(1-estcor)/nadj
    }
    
    

    }

if(method=="BW")
    {
    ##Varianzberechnung Bieler Williams 
    estp<-estimate
    resp <- as.numeric(status)  
    zet <- (resp - estp*weight)/nadj
    mzet <-mean(zet)
    varp <- (n/(n-1))*sum( (zet-mzet)^2 )
    #print("varp")
    #print(varp)
    varcor<-varp


    if(Y==0)
    {
    estcor <- 0.5/nadj 
#print(c(estcor,"estcor"))
    zet <- (resp - estcor*weight)/nadj
#print(c(zet,"zet"))
    mzet <-mean(zet)
#print(c(mzet,"mzet"))
    varcor <- (n/(n-1))*sum( (zet-mzet)^2 )
    #print("varcor1")
    #print(varcor)
    }


    if(Y==n)
    {
    estcor <- n-0.5/nadj 
    zet <- (resp - estcor*weight)/nadj
    mzet <-mean(zet)
    varcor <- (n/(n-1))*sum( (zet-mzet)^2 )
    }

}

if(method=="ADD1")
   {
   ##Varianzberechnung Bailer Portier mit add1
    estp <- (Y + 0.5)/(nadj+1)
    varp<-estp*(1-estp)/(nadj+1)
    varcor<-varp
   }


if(method=="ADD2")
   {
   ##Varianzberechnung Bailer Portier mit add1
    estp <- (Y + 1)/(nadj + 2)
    varp <- estp*(1-estp)/(nadj+2)
    varcor <- varp
   }

return(list(
 Y=Y,
 n=n,
 estimate=estimate,
 weight=weight,
 estp=estp,
 nadj=nadj,
 varp=varp,
 varcor=varcor,
 k=k))

}

