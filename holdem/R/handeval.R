handeval <-
function(num1,suit1){
    a1 = strflsh1(num1,suit1) 
    if(a1>0.5) return(8000000+a1)
    a1 = four1(num1)
    if(a1>0.5) return(7000000+a1)
    a1 = full1(num1)
    if(a1>0.5) return(6000000+a1)
    a1 = flush1(num1,suit1)
    if(a1>0.5) return(5000000+a1)
    a1 = straight1(num1)
    if(a1>0.5) return(4000000+a1)
    a1 = trip1(num1)
    if(a1>0.5) return(3000000+a1)
    a1 = twopair1(num1)
    if(a1>0.5) return(2000000+a1)
    a1 = onepair1(num1)
    if(a1>0.5) return(1000000+a1)
    a1 = nothing1(num1)
    return(a1)
} ## end of handeval

