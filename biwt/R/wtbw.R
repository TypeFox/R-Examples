wtbw <-
function(x,c1){
    ivec <- (abs(x)>c1)
    return((1-ivec)*(1-(x/c1)^2)^2)}

