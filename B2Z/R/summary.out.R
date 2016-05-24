summary.out <-
function(obj,cred){
return(c(quantile(obj,probs=c((100-cred)/(2*100),0.5,(100+cred)/(2*100))), mean(obj),sd(obj)))
}

