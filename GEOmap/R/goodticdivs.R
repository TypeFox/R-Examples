goodticdivs<-function(ddeg)
  {
    K = 0
    if(ddeg>1 & ddeg<=6){
      K = 1
      
    }
    if(ddeg>6 & ddeg<=10){
      K = 2
      
    }
    if(ddeg>10 & ddeg<=30){
      K = 5

    }                                
    if(ddeg>30 & ddeg<=50){
      K = 10

    }
    if(ddeg>50 & ddeg<=100){
      K = 15

    }
    if(ddeg>100 & ddeg<=200){
      K = 20
    }
    
    if(ddeg>200 & ddeg<=250){
      K = 30
    }
    if(ddeg>250){
      K = 60
    }

    return(K)


  }
