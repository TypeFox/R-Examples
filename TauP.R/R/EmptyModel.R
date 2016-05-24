EmptyModel <-
function()
{
  model = list()
  
model$rp=numeric(0)
model$year=numeric(0)
model$z=numeric(0)
model$vp=numeric(0)
model$vs=numeric(0)
model$rho=numeric(0)
model$qp=numeric(0)
model$qs=numeric(0)
model$conr=NaN
model$moho=NaN
model$d410=NaN
model$d520=NaN
model$d660=NaN
model$cmb=NaN
model$icb=NaN
model$dz=numeric(0)
model$dname=''

return(model)
}

