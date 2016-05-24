isr.index <-
function(hr,ht,time_l,time_u,delay)

{



u=seq(time_l,time_u,length.out=100)

dhr=splinefun(u-delay,hr(u-delay,deriv=2)^2)

dht=splinefun(u,ht(u,deriv=2)^2)

output=integrate(dht,time_l,time_u)$value/integrate(dhr,time_l-delay,time_u-delay)$value

return(output)

}
