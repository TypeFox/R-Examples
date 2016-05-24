bullwhip=function(method=c("MMSE","SMA","ES"),phi,L,p,alpha)
{ method <- match.arg(method)
  if (L==0){cat("L is at least one, the review period, ...\n")
             }
    else 
  {
 if (method=="MMSE"){r<-1+2*phi*(1-phi^L)*(1-phi^(L+1))/(1-phi)
                    }

else{ if (method=="SMA"){r<-1+2*(1-phi^p)*((L/p)^2+(L/p))
                         }      

      else{ if (method=="ES"){r<-1+(L*alpha)*(2*(1-phi)/(1-(1-alpha)*phi))+(L*alpha)^2*(1-phi)/((1-alpha)*(1-(1-alpha)*phi))
                              }  

             else {r<-"error"
                  } 
 
           }                                  
      }
 }
options(digits=5)
return(r)
}
