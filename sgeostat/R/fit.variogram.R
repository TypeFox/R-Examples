"fit.variogram" <-
  function (model = "exponential",v.object,nugget=0,sill=1000,range=1000,slope=1,...) 
{
  switch(model,
         exponential = fit.exponential(v.object,c0=nugget,ce=sill,ae=range,...),
         gaussian = fit.gaussian(v.object,c0=nugget,cg=sill,ag=range,...), 
         wave = fit.wave(v.object,c0=nugget,cw=sill,aw=range,...),
         linear = fit.linear(v.object,c0=nugget,cl=slope,...),
         spherical = fit.spherical(v.object,c0=nugget,cs=sill,as=range,...))
}
