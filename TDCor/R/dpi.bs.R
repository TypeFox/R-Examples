dpi.bs <-
function(gene_expr,ks_int,kd_int,delta_int,times,times2,time_step,noise,delay,NI)



{

mat_dpi=matrix(0,2,NI)

hr=splinefun(times2,c(rep(gene_expr[1],length.out=length(seq(-20,-1,time_step))),gene_expr))



for ( n in 1:NI)

{



# Diamond



ks=runif(4,min=ks_int[1],max=ks_int[2])

kd=runif(3,min=kd_int[1],max=kd_int[2])

   delta=runif(3,min=delta_int[1],max=delta_int[2])



  targ1 = prediction_L_norm(times,gene_expr,100*ks_int[2],ks[1],kd[1],delta[1],noise=noise)

   targ2 = prediction_L_norm(times,gene_expr,100*ks_int[2],ks[2],kd[2],delta[1],noise=noise)

targ3 = prediction_doubleL_norm(times,targ1,targ2,100*ks_int[2],ks[3],ks[4],kd[3],delta[2],delta[3],noise=noise)



ht1=splinefun(times2,c(rep(targ1[1],length.out=length(seq(-20,-1,time_step))),targ1))

ht2=splinefun(times2,c(rep(targ2[1],length.out=length(seq(-20,-1,time_step))),targ2))

ht3=splinefun(times2,c(rep(targ3[1],length.out=length(seq(-20,-1,time_step))),targ3))



mat_dpi[1,n]=dpi.index(ht1,ht2,ht3,time_l=times[1]+delay,time_u=times[length(times)],time_step=1,delay=delay)



# Co-regulation



targ3 = prediction_L_norm(times,targ1,100*ks_int[2],ks[3],kd[3],delta[2],noise=noise)

ht3=splinefun(times2,c(rep(targ3[1],length.out=length(seq(-20,-1,time_step))),targ3))

mat_dpi[2,n]=dpi.index(ht1,ht2,ht3,time_l=times[1]+delay,time_u=times[length(times)],time_step=1,delay=delay)



}



return(mat_dpi)



}
