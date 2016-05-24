
setMethod("getRiskFctBV", signature(risk = "asL4", biastype = "ANY"),
    function(risk) function(bias, var) bias^4+3*var^2+6*var*bias^2)

setMethod("getRiskFctBV", signature(risk = "asL1", biastype = "ANY"),
    function(risk) function(bias, var){
    v <- var^.5
    b <- bias
    w <- b/v
    2*v*dnorm(w)+(2*pnorm(w)-1)*b})
