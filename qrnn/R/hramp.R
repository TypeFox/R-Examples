hramp <-
function(x, lower, eps)
{   
    if(length(lower) > 1){
        mapply(hramp, x, lower, eps)
    } else{
        if (lower==-Inf){
            return(x)
        } else{
            return(ifelse(x>lower, huber(x-lower, eps), 0)+lower)
        }
    }
}
