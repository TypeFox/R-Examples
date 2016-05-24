overlap.index <-
function(net_pot,SplineList,times,time_step,mat_delay)

{

dim_net=dim(net_pot)[1]

mat_overlap=matrix(0,dim_net,dim_net)



iin=which(apply(net_pot,1,min)<0)

jin=which(apply(net_pot,2,min)<0)





for (i in iin)

{

for (j in jin)

{

if (net_pot[i,j]<0)

{

delay=mat_delay[i,j]

u=seq(times[1],times[length(times)]-delay,time_step)

overlap_ob=splinefun(u,min.element(SplineList[[i]](u),(SplineList[[j]](u))))

overlap_th=splinefun(u,min.element(SplineList[[i]](u),1-abs(SplineList[[i]](u+delay))))

mat_overlap[i,j]=integrate(overlap_ob,times[1],times[length(times)]-delay)$value/integrate(overlap_th,times[1],times[length(times)]-delay)$value

}

}

}



return(mat_overlap)

}
