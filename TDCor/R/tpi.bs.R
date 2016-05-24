tpi.bs <-
function(gene_expr,ks_int,kd_int,delta_int,times,times2,time_step,noise,delay,NI)

{



mat_tpi_part=matrix(0,3,NI)

hr=splinefun(times2,c(rep(gene_expr[1],length.out=length(seq(-20,-1,time_step))),gene_expr))



for ( n in 1:NI)

{

ks=runif(3,min=ks_int[1],max=ks_int[2])

   kd=runif(2,min=kd_int[1],max=kd_int[2])

   delta=runif(2,min=delta_int[1],max=delta_int[2])



# Feed-forward loop



y=prediction_multiL_ff_norm(times,gene_expr,ksb=100*ks_int[2],ks1=ks[1],ks21=ks[2],ks22=ks[3],kd1=kd[1],kd2=kd[2],delta1=delta[1],delta2=delta[2],noise=noise)

targ1=y$tar1

targ2=y$tar2



ht1=splinefun(times2,c(rep(targ1[1],length.out=length(seq(-20,-1,time_step))),targ1))

ht2=splinefun(times2,c(rep(targ2[1],length.out=length(seq(-20,-1,time_step))),targ2))



mat_tpi_part[1,n]=tpi.index(hr,ht1,ht2,time_l=times[1]+delay,time_u=times[length(times)],time_step=1,delay=delay)



# Cascade



   targ2 = prediction_L_norm(times,targ1,100*ks_int[2],ks[2],kd[2],delta[2],noise=noise)

ht2=splinefun(times2,c(rep(targ2[1],length.out=length(seq(-20,-1,time_step))),targ2))



mat_tpi_part[2,n]=tpi.index(hr,ht1,ht2,time_l=times[1]+delay,time_u=times[length(times)],time_step=1,delay=delay)



# Co-regulation



   targ2 = prediction_L_norm(times,gene_expr,100*ks_int[2],ks[2],kd[2],delta[1],noise=noise)

ht2=splinefun(times2,c(rep(targ2[1],length.out=length(seq(-20,-1,time_step))),targ2))



mat_tpi_part[3,n]=tpi.index(hr,ht1,ht2,time_l=times[1]+delay,time_u=times[length(times)],time_step=1,delay=delay)



}



return(mat_tpi_part)



}
