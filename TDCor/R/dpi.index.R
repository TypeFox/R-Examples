dpi.index <-
function(hr1,hr2,ht,time_l,time_u,time_step,delay)

{

u=seq(time_l,time_u,time_step)



a=abs(cor(hr1(u-delay),ht(u)))

c=abs(cor(hr2(u-delay),ht(u)))

b=abs(cor(hr1(u-delay),ht(u-delay)))

d=abs(cor(hr2(u-delay),ht(u-delay)))



return(a+b-c-d)

}
