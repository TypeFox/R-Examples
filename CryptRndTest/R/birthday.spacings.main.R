birthday.spacings.main <-
function(x,m=128,n=2^16,alpha=0.05,lambda,num.class=10){
  bds=dogumGunuAraliklari(x,m)
  
  testDGA=KSADdga(bds$ee,alpha,n,m,lambda,num.class=num.class)
  
  result=list(AD.result=testDGA$sonucAD, AD.pvalue=testDGA$ADtest[1,3], AD.statistic=testDGA$ADtest[1,1],
              KS.result=testDGA$sonucKS, KS.pvalue=testDGA$KStest[2], KS.statistic=testDGA$KStest[1],
              CS.result=testDGA$sonucKK, CS.pvalue=testDGA$KKtest[2], CS.statistic=testDGA$KKtest[1],name="BDS")
    return(result)
}
