Sc <-
function(t,accrualtime,followuptime,accrualdist){
 if(accrualdist=="uniform") return((accrualtime+followuptime-t)/accrualtime)
 if(accrualdist=="increasing") return((accrualtime+followuptime-t)^2/accrualtime^2)
 if(accrualdist=="decreasing") return((1-(followuptime-t)^2/accrualtime^2))
}

