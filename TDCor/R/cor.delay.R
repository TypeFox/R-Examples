cor.delay <-
function(delay,hr,ht,time_l,time_u,time_step,type,deriv)

{

u=seq(time_l,time_u,time_step)

if (deriv)

{

if (delay>=0)

{

k1=hr(u-delay,deriv=1)

k2=ht(u,deriv=1)

}else

{

k1=hr(u,deriv=1)

k2=ht(u+delay,deriv=1)

}

}else

{

if (delay>=0)

{

k1=hr(u-delay)

k2=ht(u)

}else

{

k1=hr(u)

k2=ht(u+delay)

}



}





if (type>0)

{

output=cor(k1,k2)

}else

{

output=cor(-k1,k2)

}

return(output)

}
