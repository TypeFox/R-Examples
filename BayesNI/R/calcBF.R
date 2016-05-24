calcBF <-
function(x1,x2,n1,n2,m,w1,w2,gfun,pH0){   
 
  
  postNull<-postNullFun(x1,x2,n1,n2,m=m,w1=w1,w2=w2,gfun=gfun)
                        
  BFX1X2<-((1-postNull)*pH0)/(postNull*(1-pH0))
                        
  if(1/BFX1X2==Inf){LogBFX1X2<- -Inf}else{LogBFX1X2<-log(BFX1X2)}
  return(LogBFX1X2)  
                 }

