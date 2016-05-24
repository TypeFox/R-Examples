orscoreci <-
function(x1,n1,x2,n2,conf.level){
  px = x1/n1
  py = x2/n2
  if(((x1==0) && (x2==0)) || ((x1==n1) && (x2==n2))){
      ul = 1/0
      ll = 0
      }
  else if((x1==0) || (x2==n2)){
       ll = 0
       theta = 0.01/n2
       ul = limit(x1,n1,x2,n2,conf.level,theta,1)
     }
  else if((x1==n1) || (x2==0)){
       ul = 1/0
       theta = 100*n1
       ll = limit(x1,n1,x2,n2,conf.level,theta,0)
     }
  else{
      theta = px/(1-px)/(py/(1-py))/1.1
      ll = limit(x1,n1,x2,n2,conf.level,theta,0)
      theta=px/(1-px)/(py/(1-py))*1.1
      ul = limit(x1,n1,x2,n2,conf.level,theta,1)
   }
  cint <- c(ll, ul)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint)
  class(rval) <- "htest"
  return(rval)
}

