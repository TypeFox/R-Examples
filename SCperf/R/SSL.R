SSL=function(method=c("MMSE","SMA","ES"),phi,L,p,alpha,SL)
{  z<- qnorm(SL, mean = 0, sd = 1)
 if (method=="MMSE"){ #Calculating the variance during the LT
                      arma1<-ARMAtoMA(ar=phi, ma=0,L);
                      arma2<-c(1,arma1);

                      totalLT <- 0
                      for (i in 1:L)
                     { valsumLT<-  (sum(arma2[1:i]))^2;
                       totalLT<- totalLT + valsumLT;
                     }

                     VarLT<-totalLT;

                     return(z*(VarLT^0.5))
                    }

else{ if (method=="SMA"){var_SMA<-(1/(1-phi^2))*(L*(1+phi)/(1-phi)-2*phi*(1-phi^L)/(1-phi)^2+((L/p)^2)*(p*(1+phi)/(1-phi)-2*phi*(1-phi^p)/(1-phi)^2)-2*(L/p)*(phi*(1-phi^L)*(1-phi^p)/(1-phi)^2))
                         return(z*(var_SMA^0.5))
                         }      

      else{ if (method=="ES"){var_ES<-(1/(1-phi^2))*(L*(1+phi)/(1-phi)-2*phi*(1-phi^L)/(1-phi)^2+((alpha*L)^2)*(1+(1-alpha)*phi)/(alpha*(2-alpha)*(1-(1-alpha)*phi))-2*alpha*L*(phi-phi^(L+1))/((1-phi)*(1-(1-alpha)*phi)))

                 return(z*(var_ES^0.5))
                              }  

             else {return("error")
                  } 

           }                                  
      }                 
}
