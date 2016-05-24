estimate.delay.cor <-
function(ht,hr,inter_time,times,time_step,thr_cor,tol,delaymax,delayin)

{

delay_p=c(0,0)

cor_p=c(0,0)



type=sign(cor(hr(inter_time),ht(inter_time)))# Sign of correlation used for guessing sign of potential interaction

u1=sapply(delayin,ssq.delay,hr=hr,ht=ht,time_l=times[1],time_u=times[length(times)],time_step=time_step,type=type,deriv=FALSE)

v1=sapply(delayin,cor.delay,hr=hr,ht=ht,time_l=times[1],time_u=times[length(times)],time_step=time_step,type=type,deriv=FALSE)

u2=sapply(delayin,ssq.delay,hr=hr,ht=ht,time_l=times[1],time_u=times[length(times)],time_step=time_step,type=type,deriv=TRUE)

v2=sapply(delayin,cor.delay,hr=hr,ht=ht,time_l=times[1],time_u=times[length(times)],time_step=time_step,type=type,deriv=TRUE)





if (min(v1)<= -min(thr_cor)) 

{

w1=sapply(delayin,ssq.delay,hr=hr,ht=ht,time_l=times[1],time_u=times[length(times)],time_step=time_step,type=-type,deriv=FALSE)

w2=sapply(delayin,ssq.delay,hr=hr,ht=ht,time_l=times[1],time_u=times[length(times)],time_step=time_step,type=-type,deriv=TRUE)

}else

{

w1=rep(max(u1),length(delayin))

w2=rep(max(u2),length(delayin))

}



if (max(abs(v1))>=min(thr_cor))

{

delay_est=c()

score_f=c()



if (max(v1)>=min(thr_cor))

{

x1=local.min(u1); wmin_u1=which.min(u1)

if (length(x1)>0)

{

x1=delayin[x1]

}else

{

x1=delayin[wmin_u1]

}





x2=local.max(v1); wmax_v1=which.max(v1) 

if (length(x2)>0)

{

x2=delayin[x2]

}else

{

x2=delayin[wmax_v1]

}



x3=local.min(u2); wmin_u2=which.min(u2) 

if (length(x3)>0)

{

x3=delayin[x3]

}else

{

x3=delayin[wmin_u2]

}



x4=local.max(v2);  wmax_v2=which.max(v2) 

if (length(x4)>0)

{

x4=delayin[x4]

}else 

{

x4=delayin[wmax_v2]

}



maxminpos=delayin[c(wmin_u1,wmax_v1,wmin_u2,wmax_v2)]

np=abs(c(min(sign(min(maxminpos)),0),max(sign(max(maxminpos)),0)))



groups=make.groups(x1,x2,x3,x4,delaymax,np)

score=rep(4,length(groups))

delay_est_vector=rep(0,length(groups))



if (length(groups)>0)

{

# the following spline functions are used to score the various delay estimations

h_ssq=splinefun(delayin,(u1-min(u1,w1))/(max(u1,w1)-min(u1,w1)))

h_cor=splinefun(delayin,max(abs(v1))-v1)

h_ssq_der=splinefun(delayin,(u2-min(u2,w2))/(max(u2,w2)-min(u2,w2)))

h_cor_der=splinefun(delayin,max(abs(v2))-v2)



for (gr in 1:length(groups))

{

candidates=groups[[gr]]

if (min(candidates)*max(candidates)>=0) 

{

delay_est_vector[gr]=mean(candidates)

part_score=rep(0,4)

part_score[1]=h_ssq(delay_est_vector[gr])

part_score[2]=h_cor(delay_est_vector[gr])

part_score[3]=h_ssq_der(delay_est_vector[gr])

part_score[4]=h_cor_der(delay_est_vector[gr])

score[gr]=mean(part_score)



}else

{

if (length(candidates)==4)

{

candidates_new=candidates[substract(1:4,which.max(abs(mean(candidates)-candidates)))]

if (mean(candidates)!=median(candidates) & min(candidates_new)*max(candidates_new)>=0 & sum(sign(candidates))!=0)

{

delay_est_vector[gr]=mean(candidates_new)

part_score=rep(0,4)

part_score[1]=h_ssq(delay_est_vector[gr])

part_score[2]=h_cor(delay_est_vector[gr])

part_score[3]=h_ssq_der(delay_est_vector[gr])

part_score[4]=h_cor_der(delay_est_vector[gr])

score[gr]=mean(part_score)

}

}

}



}

}



delay_est=c(delay_est,delay_est_vector[score<tol])

score_f=c(score_f,score[score<tol])

}



if (min(v1)<= -min(thr_cor)) # if our guess concerning the sign of the interaction may have been wrong

{

x1=local.min(w1); wmin_w1=which.min(w1)

if (length(x1)>0)

{

x1=delayin[x1]

}else

{

x1=delayin[wmin_w1]

}



v1=-v1

x2=local.max(v1); wmax_v1=which.max(v1) 

if (length(x2)>0)

{

x2=delayin[x2]

}else

{

x2=delayin[wmax_v1]

}



x3=local.min(w2); wmin_w2=which.min(w2) 

if (length(x3)>0)

{

x3=delayin[x3]

}else

{

x3=delayin[wmin_w2]

}



v2=-v2

x4=local.max(v2);  wmax_v2=which.max(v2) 

if (length(x4)>0)

{

x4=delayin[x4]

}else 

{

x4=delayin[wmax_v2]

}



maxminpos=delayin[c(wmin_w1,wmax_v1,wmin_w2,wmax_v2)]

np=abs(c(min(sign(min(maxminpos)),0),max(sign(max(maxminpos)),0)))

groups=make.groups(x1,x2,x3,x4,delaymax,np)

score=rep(4,length(groups))

delay_est_vector=rep(0,length(groups))



if (length(groups)>0)

{

# the following spline functions are used to score the various delay estimations

h_ssq=splinefun(delayin,(w1-min(u1,w1))/(max(u1,w1)-min(u1,w1)))

h_cor=splinefun(delayin,max(abs(v1))-v1)

h_ssq_der=splinefun(delayin,(w2-min(u2,w2))/(max(u2,w2)-min(u2,w2)))

h_cor_der=splinefun(delayin,max(abs(v2))-v2)



for (gr in 1:length(groups))

{

candidates=groups[[gr]]

if (min(candidates)*max(candidates)>=0) 

{

delay_est_vector[gr]=mean(candidates)

part_score=rep(0,4)

part_score[1]=h_ssq(delay_est_vector[gr])

part_score[2]=h_cor(delay_est_vector[gr])

part_score[3]=h_ssq_der(delay_est_vector[gr])

part_score[4]=h_cor_der(delay_est_vector[gr])

score[gr]=mean(part_score)



}else

{

if (length(candidates)==4)

{

candidates_new=candidates[substract(1:4,which.max(abs(mean(candidates)-candidates)))]

if (mean(candidates)!=median(candidates) & min(candidates_new)*max(candidates_new)>=0 & sum(sign(candidates))!=0)

{

delay_est_vector[gr]=mean(candidates_new)

part_score=rep(0,4)

part_score[1]=h_ssq(delay_est_vector[gr])

part_score[2]=h_cor(delay_est_vector[gr])

part_score[3]=h_ssq_der(delay_est_vector[gr])

part_score[4]=h_cor_der(delay_est_vector[gr])

score[gr]=mean(part_score)

}

}

}



}

}



delay_est=c(delay_est,delay_est_vector[score<tol])

score_f=c(score_f,score[score<tol])

}



if (length(score_f)>0)

{

delay_est_ji=delay_est[delay_est>0]

score_ji=score_f[delay_est>0]

delay_est_ij=delay_est[delay_est<0]

score_ij=score_f[delay_est<0]



# Interaction j->i



if (length(delay_est_ji)>0)

{

delay_p[1]=delay_est_ji[which.min(score_ji)]

cor_p[1]=round(cor(hr(inter_time-delay_p[1]/2),ht(inter_time+delay_p[1]/2)),3)

}



# Interaction i->j



if (length(delay_est_ij)>0)

{

delay_p[2]= -delay_est_ij[which.min(score_ij)]

cor_p[2]=round(cor(ht(inter_time-delay_p[2]/2),hr(inter_time+delay_p[2]/2)),3)

}

}

}



return(c(cor_p,delay_p))



}
