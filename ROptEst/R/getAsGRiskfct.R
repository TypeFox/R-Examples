
setMethod("get.asGRisk.fct", signature(Risk = "asL1"),
          function(Risk){return( function(r,s,b){
             rb <- r*b; w <- rb/s;
             rb*(2*pnorm(w)-1)+2*s*dnorm(w)})})

setMethod("get.asGRisk.fct", signature(Risk = "asMSE"),
          function(Risk){ return(function(r,s,b){
             rb <- r*b; rb^2+s^2})})

setMethod("get.asGRisk.fct", signature(Risk = "asL4"),
          function(Risk){ return(function(r,s,b){
             rb <- r*b; rb^4+6*s^2*rb^2+3*s^4})})

