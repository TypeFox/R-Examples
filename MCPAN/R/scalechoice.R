`scalechoice` <-
function(shape=3, pm=0.2, t=1)       ####t=max time der studie
{
scale <- ((t^shape)/(-log(1-pm)))^(1/shape)
return(scale)
}

