setMethod("getRiskFctBV", signature(risk = "asGRisk", biastype = "ANY"),
    function(risk) function(bias, var)stop("not yet implemented"))

setMethod("getRiskFctBV", signature(risk = "asMSE", biastype = "ANY"),
    function(risk) function(bias, var) bias^2+var)

setMethod("getRiskFctBV", signature(risk = "asSemivar", biastype = "onesidedBias"),
    function(risk, biastype=biastype(risk))
        function(bias,var){
            v <- var^.5
            b <- if(sign(biastype)>0) bias else -bias
            w <- b/v
            return((v^2+b^2)*pnorm(w)-b*v*dnorm(w))})

setMethod("getRiskFctBV", signature(risk = "asSemivar", biastype = "asymmetricBias"),
    function(risk, biastype=biastype(risk)){
        nu1 <- nu(biastype)[1]
        nu2 <- nu(biastype)[2]
        function(bias,var){
            v <- var^.5
            b <- if(sign(biastype)>0) bias else -bias
            w <- b/v
            return((v^2+b^2)*(nu1*pnorm(w)+nu2*pnorm(-w))+(nu2-nu1)*b*v*dnorm(w))}})

