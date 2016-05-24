# Code in this file (munoz.R) kindly provided by Francois Munoz

`GST.k` <-
function (D, exact = TRUE)
{
   GST=c();
   distrib=c();
   nk=c();
   for(i in 1:length(D[,1])) distrib=append(distrib,sum(D[i,]));

   for(k in 1:length(D[1,]))
   {
     # intra-plot similarity Fintra(k)
     if(exact) Fintrak=sum((D[,k]/sum(D[,k]))*((D[,k]-1)/(sum(D[,k])-1))) else Fintrak=sum((D[,k]/sum(D[,k]))^2);

     # conditional global similarity Fglobal(k)
     if(exact) Fglobalk=sum((D[,k]/sum(D[,k]))*((distrib-1)/(sum(distrib)-1))) else Fglobalk=sum((D[,k]/sum(D[,k]))*(distrib/sum(distrib)));

     GST=c(GST,(Fintrak-Fglobalk)/(1-Fglobalk));
     nk=c(nk,sum(D[,k]));
   };

   ret=c();
   ret$GST=GST;ret$distrib=distrib;ret$nk=nk;
   return(ret);
}

`I.k` <-
function (D, exact = TRUE)
{
   res=GST.k(D,exact);
   GST=res$GST; distrib=res$distrib; nk=res$nk;
   Ne=sum(distrib);
   
   I=(Ne-nk)*(1/GST-1)/(Ne-1);
   m=I/(I+nk-1);

   ret=c();
   ret$I=I;ret$distrib=distrib;ret$m=m;
   return(ret);
}

`etienne` <-
function (theta, m, D, log.kda = NULL, give.like = TRUE, give.log = TRUE) 
{
    D <- sort(as.table(D[D>0]), decreasing = TRUE);
    J <- sum(D)
    S <- length(D)
    if (is.null(log.kda)) {
        log.kda <- logkda.polyn(D)
    }
    A <- S:J
    if (m != 1) {
        I <- m * (J - 1)/(1 - m)
        correction.factor <- sum(exp(log.kda + lgamma(theta + 
            J) - lgamma(theta + A) + lgamma(I) - lgamma(I + J) + 
            A * log(I)))
    }
    else {
        correction.factor <- 1
    }
     if (give.like) {
        if (give.log) {
            return(S * log(theta) + lgamma(theta) - lgamma(theta + J) + log(correction.factor))
        }
        else {
            return(exp(S * log(theta) + lgamma(theta) - lgamma(theta + J) + log(correction.factor)))
        }
    }
    else {
        jj <- theta.prob(theta = theta, x = D) * correction.factor
        if (give.log) {
            return(log(jj))
        }
        else {
            return(jj)
        }
    }            
}

`logkda.polyn` <-
function (D) 
{
  a<-as.vector(D)
  a<-a[a>0]
  species=length(a)
  a<-sort(as.table(a), decreasing = TRUE);
  
  for(k in 1:species)
  {
     coeff=c(0);
     if(a[k]>1) for(i in 1:(a[k]-1)) coeff=c(coeff,exp(logS1vect[i+(a[k]-1)*(a[k]-2)/2,1]+lfactorial(i-1)-lfactorial(a[k]-1)));
     coeff=c(coeff,1);
     if(k==1) poly.kda=as.polynomial(coeff) else poly.kda=prod(poly.kda,as.polynomial(coeff));
  }
 
  return(log(coef(poly.kda)[coef(poly.kda)!=0]))
}

`optimal.params.gst` <-
function (D, exact = TRUE, ci = FALSE, cint=c(0.025,0.975), nbres=100)
{
    res=I.k(D,exact);
    I=res$I;m=res$m;
    
    #Optional computation of the confidence interval, along with the median value
    if(ci)
    {
      #Recipient of the bootstrap m estimates
      msampleboot=c();
      #Recipient of the bootstrap I estimates
      Isampleboot=c();
      #Recipient of the m confidence intervals
      msampleci=array(dim=c(length(D[1,]),3));
      #Recipient of the I confidence intervals
      Isampleci=array(dim=c(length(D[1,]),3));

      #A bootstrap resampling estimation is performed on the local sample
      for(boot in 1:nbres)
      {
        Dsample=data.frame(array(data=0,dim=c(length(D[,1]),length(D[1,]))));
        rownames(Dsample)=rownames(D);
        for (k in 1:length(D[1,]))
        {
          tab=table(sample(rownames(D),prob=D[,k], size=sum(D[,k]), replace=TRUE));
          Dsample[rownames(tab),k]=tab;
        };
        res=I.k(Dsample,exact);
        Isampleboot=rbind(Isampleboot,res$I);
        msampleboot=rbind(msampleboot,res$m);
      };
      for (k in 1:length(D[1,]))
      {
        buf=quantile(Isampleboot[,k],prob=c(cint[1],0.5,cint[2]));Isampleci[k,]=as.vector(buf[1:3]);
        buf=quantile(msampleboot[,k],prob=c(cint[1],0.5,cint[2]));msampleci[k,]=as.vector(buf[1:3]);
      };
    }

    ret=c();
    ret$I=I; ret$m=m;
    if(ci)
    {
      msampleci=as.data.frame(msampleci);colnames(msampleci)=c(paste(cint[1]*100,"%",sep=""),"Median",paste(cint[2]*100,"%",sep=""));ret$mci=msampleci;
      Isampleci=as.data.frame(Isampleci);colnames(Isampleci)=c(paste(cint[1]*100,"%",sep=""),"Median",paste(cint[2]*100,"%",sep=""));ret$Ici=Isampleci;
      ret$mboot=msampleboot;
      ret$Iboot=Isampleboot;
    };
   return(ret);
}

`optimal.params.sloss` <-
function (D, nbres=100, ci = FALSE, cint=c(0.025,0.975)) 
{
    #First stage estimation of theta
    resampl=array(dim=c(nbres,length(D[1,])));
    #Resampling using function sample with prob option:
    for(k in 1:length(D[1,])) resampl[,k-1]=sample(rownames(D),size=nbres,replace=TRUE,prob=D[,k]);
    thetaresampl=array(dim=nbres);
    #Function for Ewens' estimation
    ewens <- function(theta,nbsp,n) {(nbsp-sum(theta/(theta+(1:n)-1)))^2}
    for(k in 1:nbres) thetaresampl[k]=optimize(f=ewens,interval=c(0,100000),nbsp=length(table(resampl[k,])),n=sum(table(resampl[k,])))$minimum;
    theta=mean(thetaresampl);
    #Optional computation of the confidence interval, along with the median value
    if(ci) thetaci=quantile(thetaresampl,prob=c(cint[1],0.5,cint[2]));

    #Second-stage estimation of m
    msample=array(dim=length(D[1,]));
    Isample=array(dim=length(D[1,]));
    if(ci)
    {
      #Recipient of the bootstrap m estimates
      msampleboot=array(dim=c(length(D[1,]),nbres));
      #Recipient of the bootstrap I estimates
      Isampleboot=array(dim=c(length(D[1,]),nbres));
      #Recipient of the m confidence intervals
      msampleci=array(dim=c(length(D[1,]),3));
      #Recipient of the I confidence intervals
      Isampleci=array(dim=c(length(D[1,]),3));
    };
    for (k in 1:length(D[1,]))
    {
      Dsample<-D[,k];
      log.kda<-logkda.polyn(Dsample);
      
      f <- function(m) {
            jj <- -etienne(theta = theta, m = m, D = Dsample, log.kda = log.kda)};
      res <- optimize(f, interval=c(0,1));
      msample[k]=res$minimum;
      Isample[k]=msample[k]*(sum(Dsample)-1)/(1-msample[k]); 
      
      #If ci=T, a bootstrap resampling estimation is performed on the local sample, using the former mean theta estimate 
      if(ci)
      {
        for(boot in 1:nbres)
        {
          Dsample=table(sample(rownames(D),prob=D[,k], size=sum(D[,k]), replace=TRUE));
          log.kda<-logkda.polyn(Dsample);
          res <- optimize(f, interval=c(0,1));
          msampleboot[k,boot]=res$minimum; 
          Isampleboot[k,boot]=msampleboot[k,boot]*(sum(Dsample)-1)/(1-msampleboot[k,boot]);
        };
        buf=quantile(msampleboot[k,],prob=c(cint[1],0.5,cint[2]));msampleci[k,]=as.vector(buf[1:3]);
        buf=quantile(Isampleboot[k,],prob=c(cint[1],0.5,cint[2]));Isampleci[k,]=as.vector(buf[1:3]);
      }
    }

    ret=c();
    ret$theta=theta;ret$m=msample;ret$I=Isample;
    if(ci)
    {
      ret$thetaci=thetaci;
      msampleci=as.data.frame(msampleci);colnames(msampleci)=c(paste(cint[1]*100,"%",sep=""),"Median",paste(cint[2]*100,"%",sep=""));ret$mci=msampleci;
      Isampleci=as.data.frame(Isampleci);colnames(Isampleci)=c(paste(cint[1]*100,"%",sep=""),"Median",paste(cint[2]*100,"%",sep=""));ret$Ici=Isampleci;
      ret$thetasamp=thetaresampl;ret$mboot=msampleboot;ret$Iboot=Isampleboot;
    };
    return(ret);
}


