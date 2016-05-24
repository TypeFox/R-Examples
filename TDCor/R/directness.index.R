directness.index <-
function(hr,ht,time_l,time_u,delay,...)

{ 

u=seq(time_l,time_u,length.out=300)



Lst=list(...)

if (!is.null(Lst$type))

{ type=Lst$type}else

{ type=cor(hr(u-delay),ht(u))}



if (type>0)

{

h=splinefun(u,abs(ht(u)-hr(u)))

}else

{

h=splinefun(u,abs(1-ht(u)-hr(u)))

}



h_th=splinefun(u,abs(hr(u-delay)-hr(u)))



S=integrate(h,time_l+delay,time_u)$value

S_th=integrate(h_th,time_l+delay,time_u)$value



return(S/S_th)

}
