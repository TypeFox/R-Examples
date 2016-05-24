`ojaGradient.hyperplane` <-
function(d, x){
    return(sign(round(d%*%c(x,1),digits=8))*d[seq_along(x)])
   # return(sign(d%*%c(x,1))*d[seq_along(x)]) # This gives wrong results. Thank you, Claudia!
}

