`choropleth` <-
function(sp,v,shading=auto.shading(v),...) 
{	i = shading$cols[1+findInterval(v,shading$breaks)]
	plot(sp,col=i,...) }

