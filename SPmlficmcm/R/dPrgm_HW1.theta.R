dPrgm_HW1.theta <-
function(m,p){return(2*S1(p)*(ifelse(m==0,-S0(p),0)+
                                     ifelse(m==1,2*S0(p)-1,0)+
                                     ifelse(m==2,S0(-p),0)))
                            }
