loglikediw <-
function(x,q,beta)
{
-sum(log(ddiweibull(x,q,beta)))
}

