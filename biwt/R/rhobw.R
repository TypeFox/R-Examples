rhobw <-
function(x,c1){    
ivec <- (abs(x)>c1)
return((c1^2/6)*ivec +(1-ivec)*(x^2/2-x^4/(2*c1^2)+x^6/(6*c1^4)))}

