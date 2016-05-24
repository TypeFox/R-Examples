random.walk.tests.main=function(x,B=64,Excursion=TRUE,Expansion=TRUE,Height=TRUE,alpha=0.05){
 
  data=as.vector(x)
  
  if ((length(data)%%B)!=0){
    data=data[1:round(length(data)/B)]
  }
  result1=list(AD.result.Excursion=-1, KS.result.Excursion=-1, CS.result.Excursion=-1, AD.pvalue.Excursion=-1, 
               KS.pvalue.Excursion=-1, CS.pvalue.Excursion=-1, AD.statistic.Excursion=-1, KS.statistic.Excursion=-1,
               CS.statistic.Excursion=-1) 
  result2=list(AD.result.Height=-1, KS.result.Height=-1, CS.result.Height=-1, AD.pvalue.Height=-1, 
               KS.pvalue.Height=-1, CS.pvalue.Height=-1, AD.statistic.Height=-1, KS.statistic.Height=-1,
               CS.statistic.Height=-1)
  result3=list(AD.result.Expansion=-1, KS.result.Expansion=-1, CS.result.Expansion=-1, AD.pvalue.Expansion=-1, 
               KS.pvalue.Expansion=-1, CS.pvalue.Expansion=-1, AD.statistic.Expansion=-1, KS.statistic.Expansion=-1,
               CS.statistic.Expansion=-1)      
 
  
  if (Excursion==TRUE){
    if (sum(B==c(16, 32, 64, 128, 256))!=1){
      stop("The value of B must be in the set [16, 32, 64, 128, 256]!")
    }
    eeD=Random.walk.D(data,B)  
    test1=KSADCHRY(eeD$e,1,B,alpha)     
    result1=list(AD.result.Excursion=test1$sonucAD, KS.result.Excursion=test1$sonucKS, 
                AD.pvalue.Excursion=test1$ADtest[1,3], KS.pvalue.Excursion=test1$KStest[2], 
                AD.statistic.Excursion=test1$ADtest[1,1], KS.statistic.Excursion=test1$KStest[1],
                CS.result.Excursion=test1$sonucKK,CS.pvalue.Excursion=test1$KKtest[2],
                CS.statistic.Excursion=test1$KKtest[1])   
  }
  if (Height==TRUE){
    if (sum(B==c(64, 128, 256, 512, 1024))!=1){
      stop("The value of B must be in the set [64, 128, 256, 512, 1024]!")
    }
    eeY=Random.walk.Y(data,B)      
    test2=KSADCHRY(eeY$e,2,B,alpha)
    result2=list(AD.result.Height=test2$sonucAD, KS.result.Height=test2$sonucKS,
                 AD.pvalue.Height=test2$ADtest[1,3], KS.pvalue.Height=test2$KStest[2],
                 AD.statistic.Height=test2$ADtest[1,1], KS.statistic.Height=test2$KStest[1],
                 CS.result.Height=test2$sonucKK,CS.pvalue.Height=test2$KKtest[2],
                 CS.statistic.Height=test2$KKtest[1])    
  }
  if (Expansion==TRUE){
    if (sum(B==c(32, 64, 128))!=1){
      stop("The value of B must be in the set [32, 64, 128]!")
    }
    eeG=Random.walk.G(data,B)
    test3=KSADCHRY(eeG$e,3,B,alpha)
    result3=list(AD.result.Expansion=test3$sonucAD,KS.result.Expansion=test3$sonucKS,              
                 AD.pvalue.Expansion=test3$ADtest[1,3],KS.pvalue.Expansion=test3$KStest[2],              
                 AD.statistic.Expansion=test3$ADtest[1,1],KS.statistic.Expansion=test3$KStest[1],
                 CS.result.Expansion=test3$sonucKK,CS.pvalue.Expansion=test3$KKtest[2],
                 CS.statistic.Expansion=test3$KKtest[1])    
  }
  
  result=c(result1,result2,result3,list(name="RW",Exc=Excursion,Exp=Expansion,Hei=Height))
  return(result)
  
}