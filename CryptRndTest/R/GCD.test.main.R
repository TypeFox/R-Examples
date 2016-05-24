GCD.test.main=function(x,B=32,KS=TRUE,CSQ=TRUE,AD=TRUE,JB=TRUE,test.k=TRUE,test.g=TRUE,mu,sd,alpha=0.05){
    if (sum(x==0)>0){
        stop("Input includes invalid value: 0.")
    }
    num=as.matrix(x)
    N=nrow(num)
    y=0
    z=0
    y2=0
    oklit=0
    sig.value.k=array(NA,dim=4)
    sig.value.g=array(NA,dim=2)
    
    if (B<=64){
        for (i in 1:N){
            oklit=GCD.q(num[i,1],num[i,2])
            y[i]=oklit$k
            if (oklit$g<3) {
                y2[i]=3
            }else if (oklit$g>35){
                y2[i]=35
            }else{
                y2[i]=as.numeric(oklit$g)
            }
        }
    }else {
        for (i in 1:N){
            y[i]=GCD(num[i,1],num[i,2])$k  #uses recursive algo
            #tests for g will not be applied for 128 bit numbers         
        }
        
    }
  if (test.k==TRUE){
    teorik.Normal=round(rnorm(N,mu,sd))
    if (KS==TRUE){
      sig.value.k[1]=ks.test(y,teorik.Normal,alternative = "two.sided")$p.value
      if (sig.value.k[1]<alpha){
        KS.result.k=0
      }else {
        KS.result.k=1
      }
    }
    if (CSQ==TRUE){
      nbin=min(length(tabulate(teorik.Normal)),length(tabulate(y)))
      gozl=tabulate(y,nbins=nbin)
      bekl=tabulate(teorik.Normal,nbins=nbin)
      gozl[which(gozl==0)]=10^-5
      bekl[which(bekl==0)]=10^-5
      testCSQ.p.value=pchisq(sum((gozl-bekl)^2/bekl), df=nbin-1, ncp = 0, lower.tail = FALSE, log.p = FALSE)       
      sig.value.k[2]=testCSQ.p.value
      if (sig.value.k[2]<alpha){
        CSQ.result.k=0
      }else {
        CSQ.result.k=1
      }
    }
    if (JB==TRUE){
      testJB=  jarque.bera.test(as.matrix((y-mean(y))/sd(y)))
      sig.value.k[3]=testJB$p.value
      if (sig.value.k[3]<alpha){
        JB.result.k=0
      }else {
        JB.result.k=1
      }
    }
    if (AD==TRUE){
      if (N>=10000){
        testAD=kSamples::ad.test(list(y[1:5000],teorik.Normal[1:5000]),method="asymptotic",dist=FALSE,Nsim=500)      
      }else{
        testAD=kSamples::ad.test(list(y[1:(N/2)],teorik.Normal[1:(N/2)]),method="asymptotic",dist=FALSE,Nsim=500)      
      }
      sig.value.k[4]=testAD$ad[1,3]
      if (sig.value.k[4]<alpha){
        AD.result.k=0
      }else {
        AD.result.k=1
      }
    }
  }
  
  if (test.g==TRUE){
    kuramsalDF=0
    teo=0
    teorik=0
    rasgele=runif(N,0,1)
    bsl=1
    i=0
    while ((sum(teo)<N)){
      i=i+1      
      kuramsalDF[i]=0
      j=1:i
      kuramsalDF[i]=sum(6/((pi*j)^2))
      if (kuramsalDF[i]==0){
        kuramsalDF[i]=1
      }
      if (i==1){
        teo[i]=sum(rasgele<=kuramsalDF[i])
        teorik[bsl:(bsl+teo[i]-1)]=3
        bsl=bsl+teo[i]
      }else{
        teo[i]=sum(rasgele<=kuramsalDF[i])-sum(teo[1:i-1])
        if (i>35){         
          teorik[bsl:N]=35
          teo[i]=N+1
        }else if(i<3){          
          teorik[bsl:(bsl+teo[i]-1)]=3   
          bsl=bsl+teo[i]
        }else{         
          teorik[bsl:(bsl+teo[i]-1)]=i
          bsl=bsl+teo[i]
        }       
      } 
    }
    
    if (KS==TRUE){      
      sig.value.g[1]=ks.test(y2,teorik,alternative = "two.sided")$p.value
      if (sig.value.g[1]<alpha){
        KS.result.g=0
      }else {
        KS.result.g=1
      }
    }
    if (CSQ==TRUE){
      test2=chisq.test(y2, teorik, correct = FALSE)
      sig.value.g[2]=test2$p.value
      if (sig.value.g[2]<alpha){
        CSQ.result.g=0
      }else {
        CSQ.result.g=1
      }
    }
  }
  
  if ((test.g==TRUE) & (test.k==TRUE)) {
    result=list(sig.value.k=sig.value.k,sig.value.g=sig.value.g,KS.result.k=KS.result.k,CSQ.result.k=CSQ.result.k,
                AD.result.k=AD.result.k,JB.result.k=JB.result.k,KS.result.g=KS.result.g,CSQ.result.g=CSQ.result.g,
                name="GCD",test.k=test.k,test.g=test.g)
  }else if ((test.g==TRUE) & (test.k==FALSE)) {
    result=list(sig.value.g=sig.value.g,KS.result.g=KS.result.g,CSQ.result.g=CSQ.result.g,name="GCD",test.k=test.k,test.g=test.g)
  }else if ((test.g==FALSE) & (test.k==TRUE)) {
    result=list(sig.value.k=sig.value.k,KS.result.k=KS.result.k,CSQ.result.k=CSQ.result.k,AD.result.k=AD.result.k,
                JB.result.k=JB.result.k,name="GCD",test.k=test.k,test.g=test.g)
  }
  return(result)
}