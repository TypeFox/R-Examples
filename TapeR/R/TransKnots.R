TransKnots <-
function(knots=c(seq(0,1,0.1)),ord=4, ...){
#   ***********************************************************************************************
      c(rep(min(knots),ord),knots[2:(length(knots)-1)],rep(max(knots),ord))
    }
