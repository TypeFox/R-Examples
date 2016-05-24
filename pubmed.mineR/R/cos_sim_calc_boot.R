cos_sim_calc_boot = function(data,indices){d = t(data[indices,]);cossimdata = NULL; for (i in 1:(nrow(d) - 1)) 
       {for (j in (i+1):nrow(d)) 
                     { tempa = sum(d[i,]*d[j,]) /(sqrt(sum(d[i,]^2))* sqrt(sum(d[j,]^2))) ;  cossimdata = c(cossimdata,tempa)}
       } 
 return(cossimdata)}
