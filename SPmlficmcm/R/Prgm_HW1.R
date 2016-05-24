Prgm_HW1 <-
function(m,p){return(ifelse(m==0,(1/(1+exp(p)))^2,0)+
                                     ifelse(m==1,2*exp(p)/((1+exp(p))^2),0)+
                                     ifelse(m==2,exp(2*p)/((1+exp(p))^2),0))
                            }
