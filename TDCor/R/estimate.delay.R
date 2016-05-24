estimate.delay <-
function(dataset,tar,reg,times,time_step,thr_cor,tol,delaymax,delayspan,...)
{
	Lst=list(...)
	if( !is.null(Lst$make.graph) )
    	{ make.graph <- Lst$make.graph}else
   	{ make.graph <- TRUE}

	if( !is.null(Lst$tar.name) )
    	{ tar.name <- Lst$tar.name}else
   	{ tar.name <- tar} 	

	if( !is.null(Lst$reg.name) )
    	{ reg.name <- Lst$reg.name}else
   	{ reg.name <- reg}

	if( !is.null(Lst$main) )
    	{ maint <- Lst$main}else
   	{ maint <- paste(reg.name," and ",tar.name)} 				

	tarp=norm.data(dataset[tar,])
	regp=norm.data(dataset[reg,])
	times2=c(seq(times[1]-20*time_step,times[1]-time_step,time_step),times)
	times3=seq(times[1]-20*time_step,times[length(times)],time_step)

	ht=splinefun(times2,c(rep(tarp[1],20),tarp))
	ht=splinefun(times3,norm.data(ht(times3)))
	
	hr=splinefun(times2,c(rep(regp[1],20),regp))
	hr=splinefun(times3,norm.data(hr(times3)))
	
	inter_time=seq(times[1],times[length(times)],time_step) 		# vector of time points analysed
	delayin=seq(-delayspan,delayspan,time_step)				# vector of delays to be screened
	delay_p=c(0,0)
	cor_p=c(0,0)
	groupp=list()
	groupn=list()

	type=sign(cor(hr(inter_time),ht(inter_time)))	# Sign of correlation used for guessing sign of potential interaction
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


			if (length(groups)>0 & make.graph==T)
			{
				dev.new(width=6,height=10)
				layout(c(1,2,3,4))
				par(mar= c(4, 4, 2.5, 4))
				w=seq(times[1],times[length(times)],time_step)
				if (type>0)
				{
					maint=paste(maint," (+)",sep="")
				}else
				{
					maint=paste(maint," (-)",sep="")
				}

				plot(w,hr(w),type="l",col="black",lwd=2,ylim=c(0,1),ylab=paste("norm. transcripts ",reg.name),xlab="time (h)",main=maint)
				par(new=T)
				plot(w,ht(w),type="l",col="red",lwd=2,ylim=c(0,1),ylab="",,xlab="time (h)")
				axis(side=4,col.axis="red",col="red")
				mtext(side = 4, line=2.5,paste("norm. transcripts ",tar.name),col="red",cex=0.75)

				plot(delayin,v1,type="b",lwd=2,ylim=c(-1,1),col="black",xlab="mu (h)",ylab="F1(mu)")
				par(new=T)
				plot(delayin,u1,type="b",lwd=2,ylab="",xlab="mu (h)",axes=F,col="grey")
				axis(side=4,col.axis="dark grey",col="grey")
				mtext(side = 4, line=2.5,"F2(mu)",col="dark grey",cex=0.75)
				abline(v=delayin[which.min(u1)]-0.05,col="red",lwd=2)
				abline(v=delayin[which.max(v1)]+0.05,col="purple",lwd=2)

				plot(delayin,v2,type="b",lwd=2,ylim=c(-1,1),col="black",xlab="mu (h)",ylab="")
				par(new=T)
				plot(delayin,u2,type="b",lwd=2,ylab="F3(mu)",xlab="mu (h)",col="grey",axes=F)
				axis(side=4,col.axis="dark grey",col="grey")
				mtext(side = 4, line=2.5,"F4(mu)",col="dark grey",cex=0.75)
				abline(v=delayin[which.min(u2)]-0.05,col="red",lwd=2)
				abline(v=delayin[which.max(v2)]+0.05,col="purple",lwd=2,xlab="mu")

				scx=seq(delayin[1],delayin[length(delayin)],time_step/100)
				scy=(h_ssq(scx)+h_cor(scx)+h_ssq_der(scx)+h_cor_der(scx))/4
				plot(scx,scy,type="l",lwd=2,ylim=c(0,max(0.5,max(scy))),col="dark green",xlab="mu (h)",ylab="Total score")
				sxin=scx[scy<tol]
				syin=scy[scy<tol]
				for (i in 1:length(sxin))
				{
					lines(c(sxin[i],sxin[i]),c(syin[i],tol),col="darkseagreen3")
				}
				abline(h=tol,col="darkseagreen4",lwd=2)
				par(new=T)
				plot(scx,scy,type="l",lwd=2,ylim=c(0,max(0.5,max(scy))),col="dark green",xlab="mu (h)",ylab="Total score")

				if (length(delay_est_vector)>0)
				{		
					for (i in 1:length(delay_est_vector))
					{
						if (score[i]<tol)
						{
							abline(v=delay_est_vector[i],col="dark grey",lwd=2)
						}else
						{
							abline(v=delay_est_vector[i],col="light grey",lwd=2)
						}
					}
				}
			}

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
			if (type>0)
			{groupn=groups}else
			{groupp=groups}
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

			if (length(groups)>0 & make.graph==T)
			{
				dev.new(width=6,height=10)
				layout(c(1,2,3,4))
				par(mar= c(4, 4, 2.5, 4))
				w=seq(times[1],times[length(times)],time_step)
				if (type>0)
				{
					maint=paste(maint," (-)",sep="")
				}else
				{
					maint=paste(maint," (+)",sep="")
				}

				plot(w,hr(w),type="l",col="black",lwd=2,ylim=c(0,1),ylab=paste("norm. transcripts ",reg.name),xlab="time (h)",main=maint)
				par(new=T)
				plot(w,ht(w),type="l",col="red",lwd=2,ylim=c(0,1),ylab="",,xlab="time (h)")
				axis(side=4,col.axis="red",col="red")
				mtext(side = 4, line=2.5,paste("norm. transcripts ",tar.name),col="red",cex=0.75)

				plot(delayin,v1,type="b",lwd=2,ylim=c(-1,1),col="black",xlab="mu (h)",ylab="F1(mu)")
				par(new=T)
				plot(delayin,w1,type="b",lwd=2,ylab="",xlab="mu (h)",axes=F,col="grey")
				axis(side=4,col.axis="dark grey",col="grey")
				mtext(side = 4, line=2.5,"F2(mu)",col="dark grey",cex=0.75)
				abline(v=delayin[which.min(u1)]-0.05,col="red",lwd=2)
				abline(v=delayin[which.max(w1)]+0.05,col="purple",lwd=2)

				plot(delayin,v2,type="b",lwd=2,ylim=c(-1,1),col="black",xlab="mu (h)",ylab="")
				par(new=T)
				plot(delayin,w2,type="b",lwd=2,ylab="F3(mu)",xlab="mu (h)",col="grey",axes=F)
				axis(side=4,col.axis="dark grey",col="grey")
				mtext(side = 4, line=2.5,"F4(mu)",col="dark grey",cex=0.75)
				abline(v=delayin[which.min(u2)]-0.05,col="red",lwd=2)
				abline(v=delayin[which.max(w2)]+0.05,col="purple",lwd=2,xlab="mu")

				scx=seq(delayin[1],delayin[length(delayin)],time_step/100)
				scy=(h_ssq(scx)+h_cor(scx)+h_ssq_der(scx)+h_cor_der(scx))/4
				plot(scx,scy,type="l",lwd=2,ylim=c(0,max(0.5,max(scy))),col="dark green",xlab="mu (h)",ylab="Total score")
				sxin=scx[scy<tol]
				syin=scy[scy<tol]
				for (i in 1:length(sxin))
				{
					lines(c(sxin[i],sxin[i]),c(syin[i],tol),col="darkseagreen3")
				}
				abline(h=tol,col="darkseagreen4",lwd=2)
				par(new=T)
				plot(scx,scy,type="l",lwd=2,ylim=c(0,max(0.5,max(scy))),col="dark green",xlab="mu (h)",ylab="Total score")

				if (length(delay_est_vector)>0)
				{		
					for (i in 1:length(delay_est_vector))
					{
						if (score[i]<tol)
						{
							abline(v=delay_est_vector[i],col="dark grey",lwd=2)
						}else
						{
							abline(v=delay_est_vector[i],col="light grey",lwd=2)
						}
					}
				}
			}


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
	}else
	{
		message("The maximum correlation between the two profiles is below the threshold thr_cor.")
		message("Therefore no delay estimation was performed.")
		return()
	}

return(list(delay=delay_p[cor_p!=0],correlation=cor_p[cor_p!=0]))

}
