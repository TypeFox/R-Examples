FindDiscon <-
function(model)
{

radii=NULL


zdiff=diff(model$z)
zindy=which(zdiff==0)


radii=model$rp-model$z[zindy]
return(radii)
}

