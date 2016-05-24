.CwtData <-
function(data,reps,classes,NR,nrep,nel,L,wavelet,variance){
  x<-as.matrix(wavCWT(which(classes==1 & reps==1),wavelet=wavelet,variance=variance))
  L0<-nrow(x)*ncol(x)
  cwtdata <- mat.or.vec(L0*NR,nel)
  cwtreps<-numeric(L0*NR)
  cwtclasses<-numeric(L0*NR)
  
  cont<-0
  for (CL in 1:2){
    for (R in 1:nrep[CL]){
      cont<-cont+1
      cwtclasses[((cont-1)*L0+1):(cont*L0)]<-rep(CL,L0)
      cwtreps[((cont-1)*L0+1):(cont*L0)]<-rep(R,L0)
      for (el in 1:nel){
        cwtdata[((cont-1)*L0+1):(cont*L0),el]<-abs(as.vector(wavCWT(data[which(reps==R & classes==CL),el],wavelet=wavelet,variance=variance)))
      }
    }
  }
  return(list(cwtclasses=cwtclasses,cwtdata=cwtdata,cwtreps=cwtreps,L0=L0))

}
.easyFeaCWT <-
function()
{
  wavelet<-c()
  variance<-c()
  repeat
  {
    print("Feature type: Continuous Wavelet Transform")
    print("Should the CWT be used?")
    print("Enter Y for yes and N for no.")
    input<-readline()
    if(input%in%c("Y","y","yes","T","t","TRUE"))
    {
      feaCWT<-TRUE
      break
    } else if(input%in%c("N","n","no","F","f","FALSE"))
    {
      feaCWT<-FALSE
      break
    } else
    {
      print("Sorry! Repeat the information.")
    }
  }
  
  
  possStats<-c("haar", "gaussian1","gaussian2", "morlet")
  
  if(feaCWT==0)
  {
    print("OK. This feature shall not be used.")
  } else{
    
    repeat
    {
      print("What wavelet should be applied? Choose from:")
      print(possStats)
      wavelet<-readline()
      if(wavelet%in%possStats)
      {
        break
      } else print("Sorry! Invalid wavelet!")	
    }	
    
    repeat
    {
      print("What is the variance that should be applied?")
      variance<-as.double(readline())
      if(!is.na(variance))
      {
        if(variance>0)
        {
          break
        } else print("Sorry! You should have given a positive variance!")
      } else print("Sorry! You should have given a real number.")	
    }
    
  }
  
  result<-list(feaCWT=feaCWT, wavelets=wavelet, variance=variance)
  
  class(result)<-"feaCWT"
  
  return(result)
}
.easyFeaDoubleStatWindowing <-
function()
{
  repeat
  {
    print("Feature type: Double windowing of the original data")
    print("How many times do you want to use this kind of feature?")
    n<-as.integer(readline())
    if(!is.na(n))
    {
      break
    } else print("Sorry! You should have given an integer")
  }
  
  w<-mat.or.vec(2*n,1)
  stat<-mat.or.vec(2*n,1)
  power<-mat.or.vec(2*n,1)
  abs<-mat.or.vec(2*n,1)
  log<-mat.or.vec(2*n,1)
  mintomax<-mat.or.vec(2*n,1)
  
  possStats<-c("sum", "mean","median", "var", "prod", "min", "max", "sd", "geometric", "harmonic")
  
  if(n==0)
  {
    print("OK. This feature shall not be used.")
  } else{
    print(paste("OK. Beginning to collect the parameters of ",n," applications of this feature extraction technique:",sep=""))
  }
  
  for(i in 1:n)
  {
    if(n==0) break
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What is the desired size of the window for the first windowing process?")
      w[2*(i-1)+1]<-as.integer(readline())
      if(!is.na(w[2*(i-1)+1]))
      {
        break
      } else print("Sorry! You should have given an integer")	
    }
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What is the desired size of the window for the second windowing process?")
      w[2*i]<-as.integer(readline())
      if(!is.na(w[2*i]))
      {
        break
      } else print("Sorry! You should have given an integer")	
    }
    
    #-----------------------------------------
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What statistic should be applied for the first windowing process? Choose from:")
      print(possStats)
      stat[2*(i-1)+1]<-readline()
      if(stat[2*(i-1)+1]%in%possStats)
      {
        break
      } else print("Sorry! Invalid stat!")	
    }
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What statistic should be applied for the second windowing process? Choose from:")
      print(possStats)
      stat[2*i]<-readline()
      if(stat[2*i]%in%possStats)
      {
        break
      } else print("Sorry! Invalid stat!")	
    }	
    
    #-----------------------------------------
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What is the power that should be applied to the data in the first windowing process?")
      print("(According to the formula: data^power)")
      power[2*(i-1)+1]<-as.double(readline())
      if(!is.na(power[2*(i-1)+1]))
      {
        break
      } else print("Sorry! You should have given a real number.")	
    }
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What is the power that should be applied to the data in the second windowing process?")
      print("(According to the formula: data^power)")
      power[2*i]<-as.double(readline())
      if(!is.na(power[2*i]))
      {
        break
      } else print("Sorry! You should have given a real number.")	
    }
    
    #-----------------------------------------
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to use the absolute value of the data for the first windowing process?")
      print("(According to the formula: abs(data)^power)")
      print("Enter Y for yes and N for no.")
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        abs[2*(i-1)+1]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        abs[2*(i-1)+1]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to use the absolute value of the data for the second windowing process?")
      print("(According to the formula: abs(data)^power)")
      print("Enter Y for yes and N for no.")
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        abs[2*i]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        abs[2*i]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
    
    #-----------------------------------------
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to calculate the log of the data for the first windowing process?")
      print("(According to the formula: log(abs(data)^power+1)). Look out for negative values!!!")
      print("Enter Y for yes and N for no.")
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        log[2*(i-1)+1]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        log[2*(i-1)+1]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to calculate the log of the data for the second windowing process?")
      print("(According to the formula: log(abs(data)^power+1)). Look out for negative values!!!")
      print("Enter Y for yes and N for no.")
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        log[2*i]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        log[2*i]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
    
    #-----------------------------------------
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to sort the statistics obtained and work with quantiles for the first windowing process?")
      print("Enter Y for yes and N for no.")
      
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        mintomax[2*(i-1)+1]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        mintomax[2*(i-1)+1]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to sort the statistics obtained and work with quantiles for the second windowing process?")
      print("Enter Y for yes and N for no.")
      
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        mintomax[2*i]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        mintomax[2*i]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
  }
  
  result<-list(N=n, windowSize=w, stats=stat, 
               power=power, abs=abs, log=log, mintomax=mintomax)
  
  class(result)<-"feaDoubleStatWindowing"
  
  return(result)
}
.easyFeaPCA <-
function()
{
  feaLoadingsPCA<-0 
  feaSpecLoadingsPCA<-0
  feaSignalsPCA<-0
  n<-0
  w<-c()
  stat<-c()
  power<-c()
  abs<-c()
  mintomax<-c()
  log<-c()
  
  repeat
  {
    print("Feature type: PCA")
    print("How many principal components should be used?")
    ncomps<-as.integer(readline())
    if(!is.na(ncomps))
    {
      break
    } else print("Sorry! You should have given an integer")
  }
  
  if(ncomps==0)
  {
    print("OK. This feature shall not be used.")
  } else
  {
    repeat
    {
      print("Should the loadings of the PCA be used?")
      print("Enter Y for yes and N for no.")
      
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        feaLoadingsPCA<-1
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        feaLoadingsPCA<-0
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
    
    repeat
    {
      print("Should the loadings of the PCA applied to the spectrum be used?")
      print("Enter Y for yes and N for no.")
      
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        feaSpecLoadingsPCA<-1
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        feaSpecLoadingsPCA<-0
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
    
    
    
    
    
    
    #####
    
    repeat
    {
      print("Feature type: Windowing on PCA principal components")
      print("How many times do you want to use this kind of feature?")
      n<-as.integer(readline())
      if(!is.na(n))
      {
        break
      } else print("Sorry! You should have given an integer")
    }
    
    w<-mat.or.vec(n,1)
    stat<-mat.or.vec(n,1)
    power<-mat.or.vec(n,1)
    abs<-mat.or.vec(n,1)
    log<-mat.or.vec(n,1)
    mintomax<-mat.or.vec(n,1)
    
    possStats<-c("sum", "mean","median", "var", "prod", "min", "max", "sd", "geometric", "harmonic")
    
    if(n==0)
    {
      print("OK. This feature shall not be used.")
    } else{
      print(paste("OK. Beginning to collect the parameters of ",n," applications of this feature extraction technique:",sep=""))
    }
    
    for(i in 1:n)
    {
      if(n==0) break
      repeat
      {
        print(paste("Parameters for the application ",i,sep=""))
        print("What is the desired size of the window?")
        w[i]<-as.integer(readline())
        if(!is.na(w[i]))
        {
          break
        } else print("Sorry! You should have given an integer")	
      }
      
      repeat
      {
        print(paste("Parameters for the application ",i,sep=""))
        print("What statistic should be applied? Choose from:")
        print(possStats)
        stat[i]<-readline()
        if(stat[i]%in%possStats)
        {
          break
        } else print("Sorry! Invalid stat!")	
      }	
      
      repeat
      {
        print(paste("Parameters for the application ",i,sep=""))
        print("What is the power that should be applied to the data?")
        print("(According to the formula: data^power)")
        power[i]<-as.double(readline())
        if(!is.na(power[i]))
        {
          break
        } else print("Sorry! You should have given a real number.")	
      }
      
      repeat
      {
        print(paste("Parameters for the application ",i,sep=""))
        print("Do you want to use the absolute value of the data?")
        print("(According to the formula: abs(data)^power)")
        print("Enter Y for yes and N for no.")
        input<-readline()
        
        if(input%in%c("Y","y","yes","T","t","TRUE"))
        {
          abs[i]<-TRUE
          break
        } else if(input%in%c("N","n","no","F","f","FALSE"))
        {
          abs[i]<-FALSE
          break
        } else
        {
          print("Sorry! Repeat the information.")
        }
        
      }
      
      repeat
      {
        print(paste("Parameters for the application ",i,sep=""))
        print("Do you want to calculate the log of the data?")
        print("(According to the formula: log(abs(data)^power+1)). Look out for negative values!!!")
        print("Enter Y for yes and N for no.")
        input<-readline()
        
        if(input%in%c("Y","y","yes","T","t","TRUE"))
        {
          log[i]<-TRUE
          break
        } else if(input%in%c("N","n","no","F","f","FALSE"))
        {
          log[i]<-FALSE
          break
        } else
        {
          print("Sorry! Repeat the information.")
        }
        
      }
      
      repeat
      {
        print(paste("Parameters for the application ",i,sep=""))
        print("Do you want to sort the statistics obtained and work with quantiles?")
        print("Enter Y for yes and N for no.")
        
        input<-readline()
        
        if(input%in%c("Y","y","yes","T","t","TRUE"))
        {
          mintomax[i]<-TRUE
          break
        } else if(input%in%c("N","n","no","F","f","FALSE"))
        {
          mintomax[i]<-FALSE
          break
        } else
        {
          print("Sorry! Repeat the information.")
        }
        
      }
    }
    
    
  }
  
  result<-list(ncomps=ncomps, feaLoadingsPCA=feaLoadingsPCA, 
               feaSpecLoadingsPCA=feaSpecLoadingsPCA, 
               feaSignalsPCA=n,windowSize=w, stats=stat, 
               power=power, abs=abs, log=log, mintomax=mintomax)
  
  class(result)<-"feaPCA"
  
  return(result)
}
.easyFeaSpecStatWindowing <-
function()
{
  repeat
  {
    print("Feature type: Windowing of the spectrum")
    print("How many times do you want to use this kind of feature?")
    n<-as.integer(readline())
    if(!is.na(n))
    {
      break
    } else print("Sorry! You should have given an integer")
  }
  
  w<-mat.or.vec(n,1)
  stat<-mat.or.vec(n,1)
  power<-mat.or.vec(n,1)
  abs<-mat.or.vec(n,1)
  log<-mat.or.vec(n,1)
  mintomax<-mat.or.vec(n,1)
  
  possStats<-c("sum", "mean","median", "var", "prod", "min", "max", "sd", "geometric", "harmonic")
  
  
  if(n==0)
  {
    print("OK. This feature shall not be used.")
  } else{
    print(paste("OK. Beginning to collect the parameters of ",n," applications of this feature extraction technique:",sep=""))
  }
  
  for(i in 1:n)
  {
    if(n==0) break
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What is the desired size of the window?")
      w[i]<-as.integer(readline())
      if(!is.na(w[i]))
      {
        break
      } else print("Sorry! You should have given an integer")	
    }
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What statistic should be applied? Choose from:")
      print(possStats)
      stat[i]<-readline()
      if(stat[i]%in%possStats)
      {
        break
      } else print("Sorry! Invalid stat!")	
    }	
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What is the power that should be applied to the spectrum?")
      print("(According to the formula: spectrum^power)")
      power[i]<-as.double(readline())
      if(!is.na(power[i]))
      {
        break
      } else print("Sorry! You should have given a real number.")	
    }
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to use the absolute value of the spectrum?")
      print("(According to the formula: abs(spectrum)^power)")
      print("Enter Y for yes and N for no.")
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        abs[i]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        abs[i]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to calculate the log of the spectrum?")
      print("(According to the formula: log(abs(spectrum)^power+1)). Look out for negative values!!!")
      print("Enter Y for yes and N for no.")
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        log[i]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        log[i]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to sort the statistics obtained and work with quantiles?")
      print("Enter Y for yes and N for no.")
      
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        mintomax[i]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        mintomax[i]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
  }
  
  result<-list(N=n, windowSize=w, stats=stat, 
               power=power, abs=abs, log=log, mintomax=mintomax)
  
  class(result)<-"feaSpecStatWindowing"
  
  return(result)
}
.easyFeaStatWindowing <-
function()
{
  repeat
  {
    print("Feature type: Windowing")
    print("How many times do you want to use this kind of feature?")
    n<-as.integer(readline())
    if(!is.na(n))
    {
      break
    } else print("Sorry! You should have given an integer")
  }
  
  w<-mat.or.vec(n,1)
  stat<-mat.or.vec(n,1)
  power<-mat.or.vec(n,1)
  abs<-mat.or.vec(n,1)
  log<-mat.or.vec(n,1)
  mintomax<-mat.or.vec(n,1)
  
  possStats<-c("sum", "mean","median", "var", "prod", "min", "max", "sd", "geometric", "harmonic")
  
  if(n==0)
  {
    print("OK. This feature shall not be used.")
  } else{
    print(paste("OK. Beginning to collect the parameters of ",n," applications of this feature extraction technique:",sep=""))
  }
  
  for(i in 1:n)
  {
    if(n==0) break
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What is the desired size of the window?")
      w[i]<-as.integer(readline())
      if(!is.na(w[i]))
      {
        break
      } else print("Sorry! You should have given an integer")	
    }
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What statistic should be applied? Choose from:")
      print(possStats)
      stat[i]<-readline()
      if(stat[i]%in%possStats)
      {
        break
      } else print("Sorry! Invalid stat!")	
    }	
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("What is the power that should be applied to the data?")
      print("(According to the formula: data^power)")
      power[i]<-as.double(readline())
      if(!is.na(power[i]))
      {
        break
      } else print("Sorry! You should have given a real number.")	
    }
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to use the absolute value of the data?")
      print("(According to the formula: abs(data)^power)")
      print("Enter Y for yes and N for no.")
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        abs[i]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        abs[i]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to calculate the log of the data?")
      print("(According to the formula: log(abs(data)^power+1)). Look out for negative values!!!")
      print("Enter Y for yes and N for no.")
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        log[i]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        log[i]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
    
    repeat
    {
      print(paste("Parameters for the application ",i,sep=""))
      print("Do you want to sort the statistics obtained and work with quantiles?")
      print("Enter Y for yes and N for no.")
      
      input<-readline()
      
      if(input%in%c("Y","y","yes","T","t","TRUE"))
      {
        mintomax[i]<-TRUE
        break
      } else if(input%in%c("N","n","no","F","f","FALSE"))
      {
        mintomax[i]<-FALSE
        break
      } else
      {
        print("Sorry! Repeat the information.")
      }
      
    }
  }
  
  result<-list(N=n, windowSize=w, stats=stat, 
               power=power, abs=abs, log=log, mintomax=mintomax)
  
  class(result)<-"feaStatWindowing"
  
  return(result)
}
.fastSVM1D <-
function(features, reps, classes)
{
  num0 = reps[1]
  num1 = reps[2]
  C0<-features[,1:num0]
  C1<-features[,(num0+1):(num0+num1)]
  
  if(num0<=1) stop("There is not enougth recordings to use fast=TRUE")
  if(num1<=1) stop("There is not enougth recordings to use fast=TRUE")
  
  N0<-t(apply(C0,1,sort))
  N1<-t(apply(C1,1,sort))
  
  P0<-t(apply(C0,1,function(x) sort(x,decreasing =TRUE)))
  P1<-t(apply(C1,1,function(x) sort(x,decreasing =TRUE)))
  
  if(num0>num1)
  {
    N0<-N0[,1:num1]
    P0<-P0[,1:num1]
  }
  if(num1>num0)
  {
    N1<-N1[,1:num0]
    P1<-P1[,1:num0]
  }
  
  CusumP0<-t(apply(P0,1,cumsum))
  CusumP1<-t(apply(P1,1,cumsum))
  CusumN0<-t(apply(N0,1,cumsum))
  CusumN1<-t(apply(N1,1,cumsum))
  
  Rows=nrow(CusumP0)
  Cols=ncol(CusumP0)
  CusumP0<-cbind(rep(0,Rows),CusumP0[,-Cols])
  CusumP1<-cbind(rep(0,Rows),CusumP1[,-Cols])
  CusumN0<-cbind(rep(0,Rows),CusumN0[,-Cols])
  CusumN1<-cbind(rep(0,Rows),CusumN1[,-Cols])
  
  W=2/(N0-P1)
  Winv=2/(N1-P0)
  
  vec<-0:(Cols-1)
  E=t(2*vec-t(W*(CusumN0-CusumP1)))
  Einv=t(2*vec-t(Winv*(CusumN1-CusumP0)))
  
  Boo = (E>=0)&(W>=0)
  Booinv = (Einv>=0)&(Winv>=0)
  
  L=mat.or.vec(Rows,Cols)
  Linv=mat.or.vec(Rows,Cols)
  
  L[Boo]=W[Boo]^2/2+E[Boo]
  Linv[Booinv] = Winv[Booinv]^2/2+Einv[Booinv]
  
  Wfinal <- numeric(Rows) 
  Beta <- numeric(Rows)
  InvOrNot<-numeric(Rows)
  
  for(i in 1:Rows)
  {
    flag = 0
    if(sum(Boo[i,])>0){
      LMin = min(L[i,Boo[i,]])
    }else{flag=1}
    
    flaginv = 0
    if(sum(Booinv[i,])>0){
      LMininv = min(Linv[i,Booinv[i,]])
    } else {flaginv = 1}
    
    
    #------------------------------------------------
    #Cases 
    if((flag==0) && (flaginv==0))
    {
      if(LMin<LMininv ) 
      {
        Pos=which(L[i,]==LMin)[1]
        Wfinal[i] = W[i,Pos]
        Beta[i] = 1-W[i,Pos]*N0[i,Pos]
      }else{
        Pos=which(Linv[i,]==LMininv)[1]
        Wfinal[i] = Winv[i,Pos]
        Beta[i] = 1-Winv[i,Pos]*N1[i,Pos]
        InvOrNot[i] = 1
      }
    }
    
    if((flag==0) && (flaginv==1))
    {
      Pos=which(L[i,]==LMin)[1]
      Wfinal[i] = W[i,Pos]
      Beta[i] = 1-W[i,Pos]*N0[i,Pos]
    }
    
    if((flag==1) && (flaginv==0))
    {
      Pos=which(Linv[i,]==LMininv)[1]
      Wfinal[i] = Winv[i,Pos]
      Beta[i] = 1-Winv[i,Pos]*N1[i,Pos]
      InvOrNot[i] = 1
    }
    
    if((flag==1) && (flaginv==1))
    {
      stop("Error: The algorithm cannot converge.")
    }
  }
  
  return(list(Wfinal=Wfinal,Beta=Beta,InvOrNot=InvOrNot,classes=classes))
}
.feaSelect <-
function(P,Ptest,n.rec,Alpha, AlphaCorr,n.recTest,minacc, fast=FALSE){
  
  nrep=n.rec
  nreptest=n.recTest
  
  if(ncol(P)!=sum(nrep)) stop("Invalid parameter 'n.rec'.")
  if(ncol(Ptest)!=sum(nreptest)) stop("Invalid parameter n.recTest.")
  
  nFeatures<-nrow(P)
  svmClassRate<-rep(NA,nFeatures)	
  
  print("calculating statistics")
  nA<-nrep[1]
  nB<-nrep[2]
  if(nA<2 || nB<2) stop("There should be at least two samples of each class.\nIf you are using the 'featureSelection' function you can try changing the parameter 'testProp' to avoid this problem.\nIf you are using the 'FeatureEEG' function you can try changing the parameter 'nselec' to avoid this problem.")
  D.A <- rowSums(P[,1:nA])/nA
  D.B <- rowSums(P[,(nA+1):(nA+nB)])/nB
  Var.A.cwt <- rowSums(P[,1:nA]^2)
  Var.B.cwt <- rowSums(P[,(nA+1):(nA+nB)]^2)
  Var.A.cwt <- Var.A.cwt-nA*D.A^2
  Var.B.cwt <- Var.B.cwt-nB*D.B^2
  sp<-sqrt((Var.A.cwt+Var.B.cwt)/(nA+nB-2))
  t <- sqrt((nA)*(nB)/(nA+nB))*abs(D.A-D.B)/sp    #Teste T.
  
  print("FDR")
  pvalue <- 1-pt(abs(t),nA+nB-2)
  ord <- order(pvalue)
  ordpvalue <- pvalue[ord]
  LT<-length(t)
  max<-max(c(0,which(ordpvalue<c(1:LT)/LT*Alpha)))	
  if (max==0) {whi<-c(); return(whi)} 
  whi<-ord[1:max] 
  
  LW<-length(whi)
  
  if (LW>100) {
    st <- sort(t)
    quantil<-min(tail(st,100))
    whi<-which(t>=quantil)
    LW<-length(whi)
  }
  
  
  print("SVM")
  
  labeltrue<-factor(LETTERS[c(rep(1,nrep[1]),rep(2,nrep[2]))])
  
  ####### RESAMPLING
  if(nA<nB)
  {
    Qnt= nB-nA
    resamp = c(1:nA,sample(1:nA,nB-nA,replace=TRUE),c((nA+1):(nA+nB)))
    nrep=c(nB,nB)
  }

  if(nB<nA)
  {
    Qnt= nA-nB
    resamp = c(1:nA,sample((nA+1):(nA+nB),nA-nB,replace=TRUE),c((nA+1):(nA+nB)))
    nrep=c(nA,nA)
  } 
  
  if(nB==nA) resamp = 1:(nA+nB)
  ###### RESAMPLING
  
  
  label<-factor(LETTERS[c(rep(1,nrep[1]),rep(2,nrep[2]))])
  labeltest<-factor(LETTERS[c(rep(1,nreptest[1]),rep(2,nreptest[2]))])
  NRtest<-sum(nreptest)
  
  F<-P[whi,resamp]
  Ftest<-Ptest[whi,]
  
  whi2<-c()
  if (length(whi)==1) {F<-t(as.matrix(F)) ;Ftest<-t(as.matrix(Ftest)) }
  
  
  if(fast)
  {
    
    control<-mat.or.vec(LW,sum(nreptest))
    control[,(nreptest[1]+1):sum(nreptest)]<-1
    
    models = .fastSVM1D(F, nrep, c(0,1))
    Predictions<-.predFastSVM1D(Ftest,models)
    rates = rowMeans(Predictions==control)
    svmClassRate[whi] <- rates 
    
    for(jj in 1:length(whi))
    {
      if (rates[jj]>minacc) whi2<-c(whi2,whi[jj])  
    }
    
    
  }else{
    for (jj in 1:length(whi)){
      
      model<- svm(F[jj,], label,method= "C-classification",kernel="linear",cost=1)
      c=as.data.frame(table(predict(model, Ftest[jj,]),labeltest ))
      Taxa<-(c[1,3]+c[4,3])/NRtest
      svmClassRate[whi[jj]]<-Taxa
      if (Taxa>minacc) whi2<-c(whi2,whi[jj])
    }    
  }
  
  whi<-whi2
  if (length(whi)==1) print("1 feature selected")
  if (length(whi)<=1) return(list(Selected = whi,FDRscore = 1-2*pvalue, SVMscore = svmClassRate))
 
  
  
  print("Correlation test")
  Cor<-cor(t(cbind(P[whi2,which(labeltrue=="A")],Ptest[whi2,which(labeltest=="A")])))
  nA<-length(which(labeltrue=="A"))+length(which(labeltest=="A"))
  nB<-length(which(labeltrue=="B"))+length(which(labeltest=="B"))
  Cor2<-cor(t(cbind(P[whi2,which(labeltrue=="B")],Ptest[whi2,which(labeltest=="B")])))
  Cor<-(Cor+Cor2)/2
  LW<-length(whi)
  
  
  Tcor<-abs(log((1+Cor)/(1-Cor))/2)
  quantil<-qnorm(AlphaCorr+(1-AlphaCorr)/2,mean=0.5493,sd=sqrt(1/((nA+nB)/2-3))) 
  discart<-c()
  for (ii in 1:(LW-1))
  {
    for (jj in (ii+1):LW)
    {
      
      if (Tcor[ii,jj]>quantil) 
      {
        if (t[whi[ii]]<t[whi[jj]]) discart<-c(discart,ii)  else discart<-c(discart,jj)
      }
    }
  }
  if (length(discart)==0) whi2<-whi else whi2<-whi[-unique(discart)]
  if (is.null(whi2)) whi2<-whi[which(t[whi]==max(t[whi]))]
  
  print(paste(length(whi2)," features selected",sep=""))
  
  return(list(Selected = whi2,FDRscore = 1-2*pvalue, SVMscore = svmClassRate))
  
}
.organize <-
function(data, reps, classes, nrep, NR){
  if (is.null(ncol(data))) nel<-1 else nel<-ncol(data)
  L1<- nrow(data)/NR
  P<-mat.or.vec(L1*nel,NR)
  for (CL in 1:2){
    for (R in 1:nrep[CL]){
      P[,(CL-1)*nrep[1]+R]<-as.vector(data[which(R==reps & classes==CL),])
    }
  }
  return(P)
}
.PcaData <-
function(data,reps,classes,ncomps,nrep, nel,NR,L){
  
  pcadata<- mat.or.vec(nrow(data),ncomps)
  pcareps<- reps
  pcaclasses<- classes
  loaddata<- mat.or.vec(NR*nel,ncomps)
  loadreps<- numeric(NR*nel)
  loadclasses<- numeric(NR*nel)
  
  cont<-0
  for (CL in 1:2){
    for (R in 1:nrep[CL]){
      cont<-cont+1
      loadclasses[((cont-1)*nel+1):(cont*nel)]<-rep(CL,nel)
      loadreps[((cont-1)*nel+1):(cont*nel)]<-rep(R,nel)
      
      prov <- data[which(reps==R & classes==CL),]
      rot<-eigen(var(prov))$vectors[,1:ncomps]
      loaddata[((cont-1)*nel+1):(cont*nel),]<-rot
      pcadata[((cont-1)*L+1):(cont*L),]<-as.matrix(prov)%*%rot
    }
  }

  return(list(pcadata=pcadata,loaddata=loaddata,loadclasses=loadclasses,
  loadreps=loadreps,pcaclasses=pcaclasses,pcareps=pcareps))

}
.predFastSVM1D <-
function(features,model)
{
  Wfinal = model$Wfinal
  Beta = model$Beta 
  InvOrNot = model$InvOrNot
  classes = model$classes
  Cols = ncol(features)
  Rows = nrow(features)
  Classification <- mat.or.vec(Rows ,Cols)
  
  
  for(i in 1:Cols)
  {
    for(j in 1:Rows)
    {
      Ft = features[j,i]*Wfinal[j]+Beta[j]
      if(InvOrNot[j]==0)
      {
        if(Ft>=0) 
        {
          Classification[j,i]= classes[1]
        } else
        {
          Classification[j,i]= classes[2]
        }				
      }
      if(InvOrNot[j]==1)
      {
        if(Ft>=0) 
        {
          Classification[j,i]= classes[2]
        } else
        {
          Classification[j,i]= classes[1]
        }				
      }		
      
    }
  }
  return(Classification)
}
.Random.seed <-
c(403L, 10L, -1209790701L, 1521127601L, 1984705376L, 1536876798L, 
-919121783L, 604602299L, 287508874L, -2029668244L, -1024697873L, 
867799589L, -1346451748L, -1195287374L, 1868246957L, 1575498599L, 
-1165176690L, -245060248L, 1470410091L, -2055241975L, -1595076200L, 
65730758L, -959739359L, 1911018163L, 1124982818L, -1483163116L, 
-569016681L, -415928179L, 1369965988L, 1335658378L, -726834187L, 
-1802122673L, -2056663722L, -1613928448L, 402012483L, 1254388961L, 
-1211937200L, 79418030L, -1876216231L, -1292829077L, 842658778L, 
1234445468L, 844148383L, -321493771L, 605816204L, -304152478L, 
1774424125L, -1952871753L, 1375777118L, 199063704L, 1802306555L, 
-633097831L, 1617531432L, 18490070L, -435632303L, 1619323395L, 
154148242L, 1966566180L, 566054375L, -1722332611L, -1193098060L, 
1342030618L, 1002678373L, -778403585L, 1535783302L, -811735280L, 
1997857715L, -806674351L, 1225962496L, 1303065822L, -895138007L, 
1209293531L, -1698319894L, 1675229580L, -1228407089L, 1891828357L, 
-1634805444L, -1763654638L, -1863914739L, -1134757817L, -1949263314L, 
-388201592L, 1963216523L, -919494551L, 2112072760L, 2025200870L, 
-652630911L, 2069853203L, 947151810L, 1062730868L, 207458935L, 
1344531309L, -398210044L, 1405593386L, -98569067L, -923442385L, 
-966979018L, -1520130464L, 941339811L, -802201663L, 819624880L, 
-1984258354L, -29549511L, -1627028917L, -1389725702L, -570764548L, 
-1447143617L, -1210476523L, 468445996L, 1326931970L, -1164285219L, 
1016885719L, -1832846146L, -338422024L, -726833957L, -365897927L, 
479576712L, -1079286218L, 592164465L, 17630371L, -349860366L, 
-278504124L, 1107811847L, -480289827L, 78214356L, 2098891770L, 
-399910075L, 855482271L, 2094404902L, -1629086288L, -1431833005L, 
1080070513L, 1787287712L, 1784932414L, 1238818633L, 280757755L, 
-113194550L, 3876652L, 144406319L, -1073867291L, -317540196L, 
800976114L, -1893096083L, -1909939545L, -1364898610L, 2114134568L, 
1691652651L, -34391991L, -1725453992L, -1106382842L, -648404511L, 
928484211L, -2122383390L, 991385428L, 1424104023L, 466442829L, 
-1088938524L, 288274378L, 123312693L, 765845647L, -1427537514L, 
-1616388800L, 1383174147L, -1459138015L, 478754064L, -432276754L, 
-408353127L, -607885525L, 1197613594L, 273539292L, 387931359L, 
60074677L, -1323508660L, 594231458L, 667525117L, 511368439L, 
1272113694L, 1873192152L, -23365573L, -661216807L, 1519056744L, 
1497004694L, -132987887L, -252906045L, 1099946450L, -1993355164L, 
-1661956185L, -615399171L, 1739732212L, 1307431898L, 230983461L, 
2030276159L, 988607942L, 552718928L, 1768828531L, -1750424943L, 
-1920935232L, -156599138L, -175198103L, 1727507611L, 299787690L, 
1581712844L, 1955195023L, 2124523717L, -1150588804L, 1347101906L, 
790510669L, 635579655L, 1777684718L, 34832840L, -1744550197L, 
1666160425L, -175038600L, 306922662L, -1615984959L, -1217724589L, 
-928367614L, -663189196L, 154403511L, 686628525L, -801304636L, 
-2064877590L, 1789204565L, -1245444625L, 945384950L, 2111085344L, 
-1056492317L, 1023455617L, -2096838672L, 1829337488L, 1254464034L, 
-129798904L, 480909964L, 1767749404L, 762335378L, -458006144L, 
-1791475180L, 540361512L, 1159292474L, 1776187856L, -2043366036L, 
-1476848716L, 1184332818L, 1760883168L, -1017609268L, -1180716832L, 
-2063800286L, 98409768L, 1399027212L, 1458742076L, 179196914L, 
-324493584L, -1926591932L, 711528024L, -210087110L, -481930768L, 
160898748L, 1256373524L, -1214300478L, 1323073632L, -96483460L, 
-347776048L, 1504034082L, 681281192L, 1477077900L, 1895464124L, 
-309582030L, -1938959808L, 5404852L, -1443509688L, -1607868998L, 
-1553214960L, -430399764L, -1912730284L, -284699374L, -161810016L, 
32939724L, -1270348192L, -531590206L, -1525873880L, 1664408492L, 
1153601340L, 1525702002L, -1543499984L, 555783460L, -1247175880L, 
-809795878L, 794256272L, 1885518588L, -419409004L, -1873702846L, 
-760119040L, -1616333764L, -1392212656L, -1679163230L, -1533231800L, 
760085260L, 848945884L, 1491517778L, 2088272000L, 1585278036L, 
569424872L, 778900730L, 2022514512L, -1861718228L, 414129524L, 
-723615854L, 194624608L, -866060916L, -555725984L, -2000362014L, 
-1349475416L, 1118347468L, 1447418876L, 994537906L, -262504976L, 
513695684L, -97705640L, 1327314298L, 1711112112L, 1151276476L, 
112845908L, -471652414L, 2014812192L, -492148548L, -286534064L, 
-1355395102L, -187402136L, -2117689268L, 1041827068L, -1973159950L, 
949841152L, -2034338444L, 1253238792L, -1910402118L, -1928632880L, 
1230139372L, 1956122004L, 1688975826L, 168829280L, -821179572L, 
1036613920L, -758720638L, -564238616L, -944787220L, -874506180L, 
666005362L, 1505676912L, -212501468L, -1811089608L, -305380198L, 
-569613616L, 1467599036L, 235669140L, 1610944130L, -1000798016L, 
-196264132L, 234389648L, 655689250L, -1700590584L, 880667916L, 
-1180574564L, 2061001234L, 1407153024L, 1177138708L, -823794136L, 
-1395159750L, 971748816L, 1098195052L, 1676473652L, -1528659310L, 
1344178016L, 494092108L, -764035104L, -653104222L, 1303460776L, 
-1350518388L, 1144747196L, 966877554L, 2011832560L, 987177028L, 
-372142760L, 939106874L, 697347952L, 986029116L, -690840812L, 
517477826L, 629040864L, -1827564676L, -179793200L, -1394219870L, 
-1244510168L, -1539057396L, 1229489084L, -927197902L, -1187171904L, 
1194832436L, -670709816L, 119149754L, 1276615056L, 1559536108L, 
-1945678252L, 1128221842L, -1682662624L, 13472204L, 2014072288L, 
-384459070L, -1930780632L, 190738092L, 41898172L, -274303502L, 
-807711440L, -806194012L, 1617072824L, -31348390L, -457993840L, 
1760673660L, -1054779116L, 900857666L, -809402624L, 1462600636L, 
-840692912L, -1366896862L, 1666221384L, -627251700L, 1688416348L, 
1020868818L, -506747392L, 1084287316L, -1630177816L, 1892865018L, 
144260304L, 482175276L, -953755916L, 1465136530L, -2091147168L, 
262448396L, 742822624L, 333993314L, 1374079656L, 243167820L, 
133590396L, -1288908494L, 2068683888L, -846322876L, -1750951464L, 
1702653050L, 2047958576L, -1311804740L, -368374956L, -781844414L, 
-186349152L, 299742908L, 1883895248L, 2065567074L, 689684456L, 
1158418642L, -1474131673L, -1236163727L, -1918447178L, 1218111828L, 
-1094662171L, -208621073L, 629307688L, 550326662L, -296346189L, 
502907157L, 1607070578L, -557383744L, -1673924567L, 854469787L, 
946710892L, -550305414L, 1931301551L, 900229561L, -1977945138L, 
1501899548L, -278445427L, -246641929L, 210591232L, 994310238L, 
1234809643L, -1541480403L, -999028582L, 365583448L, 116313025L, 
425272435L, 1616588180L, -1277832158L, -1817607369L, 1814734337L, 
-1607833850L, 1803475780L, 423530101L, -1501138977L, -1940118728L, 
-1306204906L, 745134371L, -1815239963L, 648661762L, -1357938320L, 
1033922841L, -34445045L, -318658052L, -1958837718L, 1672581727L, 
-110571031L, -2117928194L, 219067116L, -1233526211L, -1850457593L, 
-827554320L, -1312184946L, 1950661563L, 699733405L, -804373686L, 
-1708703896L, 1867395409L, -325279165L, 1049798468L, -1686074382L, 
704463303L, -49727599L, -337412330L, 1326590836L, -365399547L, 
1081066063L, -2002657720L, 1381834726L, -1609912045L, 95641205L, 
1028576338L, 1763122848L, 565582729L, 1116039099L, 1246324748L, 
-723589734L, 1136821519L, -382637223L, -1810680338L, 1193611388L, 
-1038633747L, -1441504361L, -1592361888L, -708575682L, -1257072885L, 
1588997517L, -1121102150L, 1632187640L, 1039155873L, 150381139L, 
1969519988L, -1100991742L, 1237529495L, -1701284511L, 609678758L, 
1364314020L, 1106027093L, 443495167L, 893589400L, 2114689206L, 
81775939L, -1175434107L, -2016103134L, 1994874384L, -149643847L, 
1985191531L, 1744564572L, 2105034890L, 886702207L, -1377062455L, 
1695412830L, -1556010740L, -1753674531L, 792938983L, 438497808L, 
-1626460370L, 2056442843L, -1668078403L, -1191537110L, 964496072L, 
-695324687L, -936058525L, 532795748L, 111973522L, -1790465689L, 
-255612111L, -1838575498L, -142194284L, -828674651L, 236349103L, 
-93757592L, -272272570L, 1908605683L, 16948309L, -1960144846L, 
199200256L, 876969705L, 620215643L, -1266143316L, -1609583686L, 
715959663L, 597477113L, 948567950L, 63052380L, -1667827507L, 
719628343L, -377277760L, -1961253474L, -1009670549L, 1769668077L, 
-675410854L, -280577256L, -962928511L, 920737331L, 920579156L, 
-1799354270L, 1360903927L, -1903357631L, 17205190L, 1692187908L, 
1411986613L, -371567713L, -1248289800L, -1995778602L, -975133716L
)
.spec.pgram <-
function (x, taper = 0.1,fast = TRUE,detrend = TRUE, na.action = na.fail)
{
    ## Estimate spectral density from (smoothed) periodogram.
    x <- na.action(as.ts(x))
    xfreq <- frequency(x)
    x <- as.matrix(x)
    N <- N0 <- nrow(x)
    nser <- ncol(x)

        t <- 1L:N - (N + 1)/2
        sumt2 <- N * (N^2 - 1)/12
        for (i in 1L:ncol(x))
            x[, i] <- x[, i] - mean(x[, i]) - sum(x[, i] * t) * t/sumt2
    
    ## apply taper:
    x <- spec.taper(x, taper)
    ## to correct for tapering: Bloomfield (1976, p. 194)
    ## Total taper is taper*2
    u2 <- (1 - (5/8)*taper*2)
    u4 <- (1 - (93/128)*taper*2)
    x <- rbind(x, matrix(0, nrow = (nextn(N)  - N), ncol = ncol(x)))
    N <- nrow(x)
    Nspec <- floor(N/2)
    xfft <- mvfft(x)
    spec <- matrix(NA, nrow = Nspec, ncol = nser)

    for (j in 1L:ncol(x)) { # N0 = #{non-0-padded}
        spec[, j] <- Re(xfft[2:(Nspec+1), j] * Conj(xfft[2:(Nspec+1), j])/(N0*xfreq))/u2
    }
    return(spec)
}
.SpecData <-
function(data,reps,classes,NR,nrep,nel,L){
  L0<-length(.spec.pgram(1:L))
  specdata <- mat.or.vec(L0*NR,nel)
  specreps<-numeric(L0*NR)
  specclasses<-numeric(L0*NR)
  cont<-0
  for (CL in 1:2){
    for (R in 1:nrep[CL]){
      cont<-cont+1
      specclasses[((cont-1)*L0+1):(cont*L0)]<-rep(CL,L0)
      specreps[((cont-1)*L0+1):(cont*L0)]<-rep(R,L0)
      specdata[((cont-1)*L0+1):(cont*L0),]<-.spec.pgram(data[which(reps==R & classes==CL),])
    }
  }
  return(list(specclasses=specclasses,specdata=specdata,specreps=specreps,L0=L0))
}
.WinData <-
function(data, reps, classes, win, stat, power, abs, log, L, nel, nrep,mintomax) {

    Tot <- length(reps)/L
    dados<-data
    if (abs) dados <- abs(dados)
    dados <- dados^power
    if (min(dados) < 0 & log) 
    {
    	warning("Parameter log: \nSome negative values were found. Using abs=T instead.")
    	dados <- abs(dados)
    }
    
    if (log) dados <- log(dados+1)
    L2 <- L - win + 1
        Reps <- numeric(L2 * Tot)
        Classes <- numeric( L2 * 
            Tot)

    if (win>1) {
      Dados <- mat.or.vec( L2 ,nel)
      Dados2 <- mat.or.vec(L2 * Tot,nel)
      cont <- 0
      for (CL in 1:2) {
        for (R in 1:nrep[CL]) {
          cont <- cont + 1
          prov <- dados[which(classes == CL & reps == R), ]
          for (i in 1:L2) {
            Reps[(cont - 1) * L2 + i] <- R
            Classes[(cont - 1) * L2 + i] <- CL
            if (stat == "sum") Dados[i,] <-   colSums(prov[i:(i +     win - 1),])
            if (stat == "mean") Dados[ i,] <-   colMeans(prov[i:(i +     win - 1),])
            if (stat == "var") Dados[ i,] <-   colVars(prov[i:(i +     win - 1),])
            if (stat == "sd") Dados[i,] <-   sqrt(colVars(prov[i:(i +     win - 1),]))
            if (stat == "max") Dados[i,] <-   colMaxs(prov[i:(i +     win - 1),])
            if (stat == "min") Dados[i,] <-   colMins(prov[i:(i +     win - 1),])
            if (stat == "median") Dados[i,] <-   sapply(1:nel, function(g) median(prov[i:(i + win - 1),g]))
            if (stat == "prod") Dados[i,] <-   sapply(1:nel, function(g) prod(prov[i:(i + win - 1),g]))
            if (stat == "geometric") Dados[i,] <-   sapply(1:nel, function(g) (prod(prov[i:(i + win - 1),g]))^(1/win))
            if (stat == "harmonic") Dados[i,] <-   sapply(1:nel, function(g) 1/mean(1/prov[i:(i + win - 1),g]))     
          }
          if (mintomax) Dados2[((cont - 1) * L2 + 1):(cont*L2),]<-sapply(1:nel, function(g) sort(Dados[,g])) else {
            Dados2[((cont - 1) * L2 + 1):(cont*L2),]<- Dados
          }#mintomax
  
        }
      }

    } else {
      Dados2 <- dados
      Reps <- reps
      Classes<-classes
      if (mintomax) {
        Dados <- mat.or.vec( L2 ,nel)
        Dados2 <- mat.or.vec(L2 * Tot,nel)
        cont <- 0
        for (CL in 1:2) {
          for (R in 1:nrep[CL]) {
              cont <- cont + 1
              Dados <- dados[which(classes == CL & reps == R), ]
          
              Dados2[((cont - 1) * L2 + 1):(cont*L2),]<-sapply(1:nel, function(g) sort(Dados[,g])) 
            }
        }

      } #if mintomax
 
    }

    return(list(windata=Dados2, winreps=Reps, winclasses=Classes,L2=L2))
}
.winPrac <-
function(data,win,stat,power,abs,log,mintomax,W,nel){

  if (is.null(ncol(data))) data<-as.matrix(data)
  
  L<-nrow(data)
  
  if(mintomax){

    if (abs) data <- abs(data)
    data <- data ^power
    if (min(data) < 0 & log) 
    {
    	warning("Parameter log: \nSome negative values were found. Using abs=T instead.")
    	data <- abs(data)
    }
    
    if (log) data <- log(data +1)
    
    L2 <- L - win + 1
    Dados <- mat.or.vec( L2 ,nel)
    if(win!=1){
      for (i in 1:L2) {
          if (stat == "sum") Dados[i,] <-   colSums(data[i:(i +     win - 1),])
          if (stat == "mean") Dados[ i,] <-   colMeans(data[i:(i +     win - 1),])
          if (stat == "var") Dados[ i,] <-   colVars(data[i:(i +     win - 1),])
          if (stat == "sd") Dados[i,] <-   sqrt(colVars(data[i:(i +     win - 1),]))
          if (stat == "max") Dados[i,] <-   colMaxs(data[i:(i +     win - 1),])
          if (stat == "min") Dados[i,] <-   colMins(data[i:(i +     win - 1),])
          if (stat == "median") Dados[i,] <-   sapply(1:nel, function(g) median(data[i:(i + win - 1),g]))
          if (stat == "prod") Dados[i,] <-   sapply(1:nel, function(g) prod(data[i:(i + win - 1),g]))
          if (stat == "geometric") Dados[i,] <-   sapply(1:nel, function(g) (prod(data[i:(i + win - 1),g]))^(1/win))
          if (stat == "harmonic") Dados[i,] <-   sapply(1:nel, function(g) 1/mean(1/data[i:(i + win - 1),g]))
         
      }
    } else {
      for (i in 1:L2) {
          if (stat == "sum") Dados[i,] <-   colSums(matrix(data[i:(i +     win - 1),],ncol=nel))
          if (stat == "mean") Dados[ i,] <-   colMeans(matrix(data[i:(i +     win - 1),],ncol=nel))
          if (stat == "var") Dados[ i,] <-   colVars(matrix(data[i:(i +     win - 1),],ncol=nel))
          if (stat == "sd") Dados[i,] <-   sqrt(colVars(matrix(data[i:(i +     win - 1),],ncol=nel)))
          if (stat == "max") Dados[i,] <-   colMaxs(matrix(data[i:(i +     win - 1),],ncol=nel))
          if (stat == "min") Dados[i,] <-   colMins(matrix(data[i:(i +     win - 1),],ncol=nel))
          if (stat == "median") Dados[i,] <-   sapply(1:nel, function(g) median(matrix(data[i:(i +     win - 1),g],ncol=nel)))
          if (stat == "prod") Dados[i,] <-   sapply(1:nel, function(g) prod(matrix(data[i:(i +     win - 1),g],ncol=nel)))
          if (stat == "geometric") Dados[i,] <-   sapply(1:nel, function(g) (prod(matrix(data[i:(i +     win - 1),g],ncol=nel)))^(1/win))
          if (stat == "harmonic") Dados[i,] <-   sapply(1:nel, function(g) 1/mean(1/matrix(data[i:(i +     win - 1),g],ncol=nel)))     
      }
    }
    Dados<-sapply(1:nel, function(g) sort(Dados[,g])) 
    Saida<-Dados[W]

  } else {
    
    LA<-length(W)
    Saida<-numeric(LA)
    
    for(ij in 1:LA) {
      vvv<-(W[ij]-1)%%(L-win+1)+1
      val<-data[vvv:(vvv+win-1),floor((W[ij]-1)/(L-win+1))+1]
       
       if (abs) val <- abs(val)
        val <- val ^power
        if (min(val) < 0 & log) 
        {
        	warning("Parameter log: \nSome negative values were found. Using abs=T instead.")
        	val <- abs(val)
        }
    
        if (log) val <- log(val +1)
    


        if (stat=="sum") Saida[ij]<-sum(val)
        if (stat=="mean") Saida[ij]<-mean(val)
        if (stat=="var") Saida[ij]<-var(val)
        if (stat=="sd") Saida[ij]<-sqrt(var(val))
        if (stat=="max") Saida[ij]<-max(val)
        if (stat=="min") Saida[ij]<-min(val)
        if (stat=="median") Saida[ij]<-median(val)
        if (stat=="prod") Saida[ij]<-prod(val)
        if (stat=="geometric") Saida[ij]<-(prod(val))^(1/win)
        if (stat=="harmonic") Saida[ij]<-1/mean(1/val)
    }

  } #else mintomax

  return(Saida)
}
.winPrac2 <-
function(data,win,stat,power,abs,log,mintomax,nel,L){

  if (is.null(ncol(data))) data<-as.matrix(data)

  if (abs) data <- abs(data)
  data <- data ^power
  if (min(data) < 0 & log) 
  {
    	warning("Parameter log: \nSome negative values were found. Using abs=T instead.")
    	data <- abs(data)
  }
    
  if (log) data <- log(data +1)
    


  L2 <- L - win + 1
  Dados <- mat.or.vec( L2 ,nel)

  for (i in 1:L2) {
        if (stat == "sum") Dados[i,] <-   colSums(data[i:(i +     win - 1),])
        if (stat == "mean") Dados[ i,] <-   colMeans(data[i:(i +     win - 1),])
        if (stat == "var") Dados[ i,] <-   colVars(data[i:(i +     win - 1),])
        if (stat == "sd") Dados[i,] <-   sqrt(colVars(data[i:(i +     win - 1),]))
        if (stat == "max") Dados[i,] <-   colMaxs(data[i:(i +     win - 1),])
        if (stat == "min") Dados[i,] <-   colMins(data[i:(i +     win - 1),])
        if (stat == "median") Dados[i,] <-   sapply(1:nel, function(g) median(data[i:(i + win - 1),g]))
        if (stat == "prod") Dados[i,] <-   sapply(1:nel, function(g) prod(data[i:(i + win - 1),g]))
        if (stat == "geometric") Dados[i,] <-   sapply(1:nel, function(g) (prod(data[i:(i + win - 1),g]))^(1/win))
        if (stat == "harmonic") Dados[i,] <-   sapply(1:nel, function(g) 1/mean(1/data[i:(i + win - 1),g]))
         
  }

  if (mintomax) Dados<-sapply(1:nel, function(g) sort(Dados[,g]))
return(Dados)
}

.feaSelectMult <-
  function(features, totVec,Alpha, AlphaCorr,minacc, fast=FALSE, testProp){
    
    nclasses = length(totVec) - 1
    P = features
    
    nFeatures<-nrow(P)
    svmClassRate<-rep(NA,nFeatures)  
    
    print("Initializing analysis of variance")
    
    
    #-------------------------------------------------------------
    #INI: Calculando as estatisticas
    N <- ncol(P)
    Y. <- mat.or.vec(nrow(P),nclasses)
    for(i in 1:nclasses)
    {
      Y.[,i] <- (rowSums(P[,(1+sum(totVec[1:i])):(sum(totVec[1:(i+1)]))]))^2
    }
    Y.. <- rowSums(P)^2/N
    SST <- rowSums(P^2)-Y..
    SSTrat <- t(apply(Y., 1, function(x) x/totVec[-1]))
    SSTrat <- rowSums(SSTrat) - Y..
    SSE <- SST - SSTrat
    
    F0 <- (N-nclasses)/(nclasses-1)*SSTrat/SSE
    #END: Calculando as estatisticas   
    #-------------------------------------------------------------
    
    #-------------------------------------------------------------
    #INI: Teste FDR
    print("FDR")
    pvalue <- 1-pf(F0,nclasses-1,N-nclasses)
    ord <- order(pvalue)
    ordpvalue <- pvalue[ord]
    LT<-length(F0)
    max<-max(c(0,which(ordpvalue<c(1:LT)/LT*Alpha)))  
    if (max==0) {whi<-c(); return(whi)} 
    whi<-ord[1:max] 
    
    LW<-length(whi)
    
    if (LW>100) {
      st <- sort(F0[whi])
      quantil<-min(tail(st,100))
      whi<-which(F0>=quantil)
      LW<-length(whi)
    }  
    #END: Teste FDR
    #-------------------------------------------------------------
    
    print("SVM")
    
    #-------------------------------------------------------------
    #INI: Calculating n.rec e n.recTest
    if(testProp<0.5)
    {  
      n.rec = floor(totVec[-1]*(1-testProp))
      n.recTest = totVec[-1] - n.rec
      n.rec = max(n.rec)
    } else
    {
      n.recTest = floor(totVec[-1]*(testProp))
      n.rec = totVec[-1] - n.recTest
      n.rec = max(n.rec)
    }
    #END: Calculating n.rec e n.recTest
    #-------------------------------------------------------------
    
    
    
    #-------------------------------------------------------------
    #INI: Calculating labeltrue do banco original
    labeltrue <- c()
    for(i in 1:nclasses)
    {
      labeltrue <- c(labeltrue,rep(i,totVec[i+1]))
    }
    labeltrue<-factor(LETTERS[labeltrue])
    #END: Calculating labeltrue do banco original
    #-------------------------------------------------------------
    
    
    
    #-------------------------------------------------------------
    #INI: Separando o banco de testes
    resamp <- c()
    for(i in 1:nclasses)
    {
      resamp<-c(resamp, sample(sum(totVec[1:i])+1:(totVec[i+1]),n.recTest[i]))
    }
    Ftest = P[whi,resamp]
    P = P[,-resamp]
    #END: Separando o banco de testes
    #-------------------------------------------------------------
    
    
    
    #-------------------------------------------------------------
    #INI: RESAMPLING
    maxRecTest <- max(n.recTest)
    maxSamp <- max(totVec)-maxRecTest
    resamp <- c()
    totVecAnt<-totVec
    totVec<-totVec-c(0,n.recTest)
    for(i in 1:nclasses)
    {
      if(totVec[i+1] != maxSamp)
      {
        resamp = c(resamp,c(sum(totVec[1:i])+ 1:totVec[i+1],sum(totVec[1:i]) + sample(1:totVec[i+1],maxSamp - totVec[i+1], replace = TRUE)))
      } else
      {
        resamp <- c(resamp,sum(totVec[1:i])+1:totVec[i+1])
      }
    }
    totVec<-totVecAnt
    P <- P[,resamp]
    #END: RESAMPLING
    #-------------------------------------------------------------
    
    
    
    #-------------------------------------------------------------
    #INI: Compondo label de treinamento e testes
    label<-c()
    labeltest<-c()
    
    for(i in 1:nclasses)
    {
      label<-c(label,rep(i,n.rec))
      labeltest<-c(labeltest,rep(i,n.recTest[i]))
    }   
    
    label<-factor(LETTERS[label])
    labeltest<-factor(LETTERS[labeltest])
    NRtest<-sum(n.recTest)
    #END: Compondo label de treinamento e testes
    #-------------------------------------------------------------
    
    
    F<-P[whi,]  
    
    whi2<-c()
    if (length(whi)==1) {F<-t(as.matrix(F)) ;Ftest<-t(as.matrix(Ftest)) }
    
    if(fast)
    {
      control<-mat.or.vec(LW,sum(n.recTest))
      control[,(n.recTest[1]+1):sum(n.recTest)]<-1
      
      models = .fastSVM1D(F, c(n.rec,n.rec), c(0,1))
      Predictions<-.predFastSVM1D(Ftest,models)
      rates = rowMeans(Predictions==control)
      svmClassRate[whi] <- rates 
      
      for(jj in 1:length(whi))
      {
        if (rates[jj]>minacc) whi2<-c(whi2,whi[jj])  
      }
      
      
    }else{
      for (jj in 1:length(whi)){
        
        model<- svm(F[jj,], label,method= "C-classification",kernel="linear",cost=1)
        c=as.data.frame(table(predict(model, Ftest[jj,]),labeltest ))
        if(nclasses==2)
        {
          Taxa<-(c[1,3]+c[4,3])/NRtest
        }else
        {
          sum = 0
          for(i in 1:nrow(c))
          {
            if(c[i,1]==c[i,2]) sum = sum + c[i,3]
          }
          Taxa<-sum/sum(c[,3])
        }
        
        svmClassRate[whi[jj]]<-Taxa
        if (Taxa>minacc) whi2<-c(whi2,whi[jj])
        
      }    
    }
    whi<-whi2
    if (length(whi)==1) print("1 feature selected")
    if (length(whi)<=1) return(list(Selected = whi,FDRscore = 1-2*pvalue, SVMscore = svmClassRate))
    
    
    
    print("Correlation test")
    LW<-length(whi)
    Cor<-mat.or.vec(LW,LW)
    for(i in 1:nclasses)
    {
      Cor2<-cor(t(features[whi2,which(labeltrue==LETTERS[i])]))
      Cor = Cor+Cor2
    }
    Cor<-Cor/nclasses
    
    
    Tcor<-abs(log((1+Cor)/(1-Cor))/2)
    quantil<-qnorm(AlphaCorr+(1-AlphaCorr)/2,mean=0.5493,sd=sqrt(1/((sum(totVec))/2-3))) 
    discart<-c()
    for (ii in 1:(LW-1))
    {
      for (jj in (ii+1):LW)
      {
        
        if (Tcor[ii,jj]>quantil) 
        {
          if (F0[whi[ii]]<F0[whi[jj]]) discart<-c(discart,ii)  else discart<-c(discart,jj)
        }
      }
    }
    if (length(discart)==0) whi2<-whi else whi2<-whi[-unique(discart)]
    if (is.null(whi2)) whi2<-whi[which(F0[whi]==max(F0[whi]))]
    
    print(paste(length(whi2)," features selected",sep=""))
    
    return(list(Selected = whi2,FDRscore = 1-2*pvalue, SVMscore = svmClassRate))
    
  }