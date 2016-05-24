psibw <-
function(x,c1){
ivec <- (abs(x)>c1)
    return((1-ivec)*(x*(1-(x/c1)^2)^2))}

