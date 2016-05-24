tpi.index <-
function(hr,ht1,ht2,time_l,time_u,time_step,delay)

{

u=seq(time_l,time_u,time_step)

a=abs(cor(ht1(u),ht2(u)))

b=abs(cor(hr(u-delay),ht2(u)))

c=abs(cor(ht1(u-delay),ht2(u)))



return(a+b-2*c)

}
