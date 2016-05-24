prateR.exact <-
function(x1,x2,e1,e2,delta.grd, midp=T)
      {
        if(x1+x2>0)
          {p=e1/(e1+delta.grd*e2)
           if(midp==T)
             {pv1=  pbinom(x1, x1+x2, p)-dbinom(x1, x1+x2, p)*0.5
              pv2=1-pbinom(x1, x1+x2, p)+dbinom(x1, x1+x2, p)*0.5
              }else{
                    pv1=  pbinom(x1, x1+x2, p)
                    pv2=1-pbinom(x1, x1+x2, p)+dbinom(x1, x1+x2, p)
                    }
           }
       if(x1+x2==0)
           pv1=pv2=rep(1, length(delta.grd))    
       return(list(pv1=pv1, pv2=pv2))
      }
