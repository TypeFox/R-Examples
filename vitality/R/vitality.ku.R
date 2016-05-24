vitality.ku<-function(time,sdata,rc.data=F,se=F,gfit=F,datatype="CUM",ttol=.000001,
                      init.params=F,lower=c(0,-1,0,0),upper=c(100,100,50,50),pplot=T,
                      tlab="days",lplot=F,cplot=F,Iplot=F,silent=F,L=0){
# 
#
#  Vitality based survival model: parameter fitting routine:    VERSION: 10/12/2007; DS 2014/11/17
#
# REQUIRED PARAMETERS: 
#    time - time component of data: time from experiment start.  Time should 
#          start after the imposition of a stressor is completed.
#    sdata - survival or mortality data.  The default expects cumulative 
#   	       survival fraction.  If providing incremental mortality fraction 
#          instead, use option: datatype="INC".
#          The default also expects the data to represent full mortality.   
#          Otherwise, use option: rc.data=T to indicate right censored data.
#
# OPTIONAL PARAMETERS:
#     rc.data =T  - specifies Right Censored data.   If the data does not 
#          represent full mortality, it is probably right censored.  The default 
#          is rc.data=F.  A third option is rc.data="TF".  Use this case to add
#          a near-term zero survival point to data which displays nearly full 
#          mortality ( <.01 survival at end).  If rc.data=F but the data does
#	      not show full mortality, rc.data="TF" will be 
#          invoked automatically. 
#     se =<population>  calculates the standard errors for the MLE parameters.  
#          Default is se=F.  The initial study population is necessary for 
#          computing these standard errors. 
#     gfit =<population>  provides a Pearson C type test for goodness of fit.  
#	     Default is gfit=F. The initial study population is necessary for 
#	     computing goodness of fit. 
#     datatype ="CUM" -cumulative survival fraction data- is the default.   
#          Other option: datatype="INC" - for incremental mortality fraction 
#          data.  ttol (stopping criteria tolerence.)  Default is .000001 . 
#          specify as ttol=.0001.
#          If one of the liklihood plots (esp. for "k") does not look optimal, 
#	     try decreasing ttol.   If the program crashes, try increasing ttol.
#     init.params =F  has the routine choose initial parameter estimates for 
#          r,s,k,u  (default: =F). If you wish to specify initial param values 
#          rather than have the routine choose them, specify 
#          init.params=c(r,s,k,u) in that order (eg. init.params=c(.1,.02,.003)).
#     pplot =T provides plots of cumulative survival and incremental mortality - 
#	     for both data and fitted curves (default: =T).  pplot=F provides no 
#	     plotting.  A third option:  pplot=n  (n>=1) extends the time axis of 
#	     the fitting plots 	(beyond the max time in data).  For example: 
#          pplot=1.2 extends the time axis by 20%.  (Note:  the incremental 
#          mortality plot is a continuous representation of the appropriately-
#          binned histogram of incremental mortalities.)
#     tlab ="<time units>" specifies units for x-axis of plots.  Default is 
#          tlab="days".
#     lplot =T provides likelihood function plotting (default =T).  
#          Note:  these plots are not "likelihood profiles" in that while one 
#          parameter is varied, the others are held fixed, rather than 
#          re-optimized. (must also have pplot=T.)
#     cplot =T provides a likelihood contour plot for a range of r and s values 
#          (can be slow so default is F).   Must also have lplot=T (and pplot=T) 
#          to get contour plots.
#     silent =T stops all print and plot options (still get most warning and all 
#          error messages) Default is F.  A third option, silent="verbose" also 
#          enables the trace setting in the ms (minimum sum) S-Plus routine. 
#     L =0 times of running simulated annealing. Default is 0, use Newton-Ralphson method only.
# RETURN:
#     vector of final MLE r,s,k,u parameter estimates.
#     standard errors of MLE parameter estimates (if se=<population> is
#     specified).
#
#  --Check/prepare Data---
	dTmp<-dataPrep(time,sdata,datatype,rc.data)
	time<-dTmp$time
	sfract<-dTmp$sfract
	x1<-dTmp$x1
	x2<-dTmp$x2
	Ni<-dTmp$Ni
	rc.data<-dTmp$rc.data

#  --Produce initial parameter values---
	if(length(init.params)==1) {  
          ii<-indexFinder(sfract,0.5)
           if (ii == -1) {
	           warning("ERROR: no survival fraction data below the .5 level.   Can not use the initial r s k u estimator.  You must supply initial r s k u estimates")
	           return(-1)
            }
            else rsk<-c(1/time[ii],0.1,0.01,0.1)  
	} else {                      # use user specified init params
		rsk <-init.params
	}
	if(rsk[1] == -1) {
		stop
	}
	if (silent == F) {
	  print(cbind(c("Initial r","initial s","initial k","initial u"),rsk))
    }

#  --create dataframe for sa---
	dtfm <- data.frame(x1=x1,x2=x2,Ni=Ni)
	#param(dtfm,"r")    
	#param(dtfm,"s")   
	#param(dtfm,"k")
	#param(dtfm,"u")
	
#  --run MLE fitting routine (simulated annealing)---   
#  --L>=1 simulated anealing will be conducted to generate the initial values--     
     if (L>=1){    
         vfit.sa.temp<-sapply(1:L,function(r0,s0,k0,u0,xx1,xx2,NNi)
            sa.lt.ku(rsk[1],rsk[2],rsk[3],rsk[4],x1,x2,Ni))

         ind<-which.min(vfit.sa.temp[5,])             
         vfit.sa<-vfit.sa.temp[,ind]

     #  --Newton-Ralphson algorithm using the results from simulated anealing as initial values--#
         fit.nlm<-nlminb(vfit.sa[1:4],objective=logLikelihood.ku,lower=lower,upper=upper,xx1=x1,xx2=x2,NNi=Ni)
     } else if (L==0){   # --conduct Newton-Ralphoson algorithm directly --
               fit.nlm<-nlminb(rsk,objective=logLikelihood.ku,lower=lower,upper=upper,xx1=x1,xx2=x2,NNi=Ni)
          } else stop("ERROR:  L should be a positive integer.")	

# --save final param estimates---
     r.final<-fit.nlm$par[1]
     s.final<-abs(fit.nlm$par[2])
     k.final<-fit.nlm$par[3]
     u.final<-fit.nlm$par[4]
     mlv<-fit.nlm$obj  
     if (silent ==F) {print(cbind(c("estimated r", "estimated s", "estimated k", "estimated u", "minimum -loglikelihood value"), c(r.final, s.final, k.final, u.final, mlv)))}
   
#  ==end MLE fitting===
   
#  --compute standard errors---
	if (se != F) {
		s.e. <-stdErr.ku(r.final, s.final, k.final, u.final, x1, x2, Ni, se)  
	      if (silent==F){print(cbind(c("sd for r","sd for s","sd for k","sd for u"), s.e.))}  
      }

# #  --compute AIC --  
#       if (AIC !=F) {
#         AIC.calc<-function(pop,value){
#           #  Routine to calculate AIC value
#           #  pop - population number of sample
#           #  value  - the minimal value of the -loglikelihood
#           return(pop*(2*4+2*value)) 
#         }
#             AIC.value<-AIC.calc(AIC,mlv)
#             if (silent==F){print(c("AIC value for model fitting:", AIC.value))}
#       }
	
#  --plotting and goodness of fit---
   if (pplot != F || gfit != F) {
	   plotting.ku(r.final,s.final,k.final,u.final,mlv,time,sfract,x1,x2,Ni,pplot,tlab,lplot,cplot,Iplot,gfit)
   }
# ............................................................................................
			
# --return final param values---
	sd<-5   #significant digits of output
	if(se != F ) {
		params<-c(r.final,s.final,k.final,u.final)
            pvalue<-c(1-pnorm(r.final/s.e.[1]),1-pnorm(s.final/s.e.[2]),1-pnorm(k.final/s.e.[3]),1-pnorm(u.final/s.e.[4]))
            std<-c(s.e.[1],s.e.[2],s.e.[3],s.e.[4])
		out<-signif(cbind(params,std,pvalue),sd)
		return(out)
	} else {
		return(signif(c(r.final,s.final,k.final,u.final),sd))
	}
}

#=SurvFn.ku================================================================================================

SurvFn.ku<-function(xx,r,s,k,u){  #  The cumulative survival distribution function.
	yy<-u^2+s^2*xx
	# pnorm is: cumulative prob for the Normal Dist.
	tmp1 <- sqrt(1/yy) * (1 - xx * r)    #  xx=0 is ok.  pnorm(+-Inf) is defined
	tmp2 <- sqrt(1/yy) * (1 + xx * r+2*u^2*r/s^2)

	# --safeguard if exponent gets too large.---
	tmp3 <- 2*r/(s*s)+2*u^2*r^2/s^4
	
	if (tmp3 >250) {   
		q <-tmp3/250 
		
		if (tmp3 >1500) {
			q <-tmp3/500
		}
		
		valueFF <-(1.-(pnorm(-tmp1) + (exp(tmp3/q) *pnorm(-tmp2)^(1/q))^(q)))*exp(-k*xx) 
	}		    
	else {
		valueFF <-(1.-(pnorm(-tmp1) + exp(tmp3) *pnorm(-tmp2)))*exp(-k*xx)   #1-G
	}
	if ( all(is.infinite(valueFF)) ) {
		warning(message="Inelegant exit caused by overflow in evaluation of survival function.  
		Check for right-censored data. Try other initial values.")
	}
	
	return(valueFF)	
}

#=survProbInc.ku===========================================================================================

survProbInc.ku<-function(r,s,k,u,xx1,xx2){ # calculates incremental survival probability
	value.iSP <--(SurvFn.ku(xx2,r,s,k,u) - SurvFn.ku(xx1,r,s,k,u))
	value.iSP[value.iSP < 1e-18] <-1e-18   # safeguards against taking Log(0)
	value.iSP
}

#=logLikelihood.ku=========================================================================================

logLikelihood.ku<-function(par,xx1,xx2,NNi){
	#returns vector of terms in log likelihood (sum them to get the log likelihood)
	# --calculate incremental survival probability--- (safeguraded >1e-18 to prevent log(0))
	iSP <-  survProbInc.ku(par[1],par[2],par[3],par[4],xx1,xx2)	
  loglklhd <--NNi*log(iSP)
	# add smooth penalty to log-likelihood if k <0
# 	if (k < 0) {
# 		loglklhd <-loglklhd + k*k*1e4
# 	}
   return(sum(loglklhd))
}

#=logLikelihood4.ku========================================================================================
#--the version of -log likelihood for function nlminb --
 
# logLikelihood4.ku<-function(par,xx1,xx2,NNi){
# 	#returns vector of terms in log likelihood (sum them to get the log likelihood)
# 	# --calculate incremental survival probability--- (safeguraded >1e-18 to prevent log(0))
# 	iSP <-  survProbInc.ku(par[1],par[2],par[3],par[4],xx1,xx2)
#   loglklhd <--NNi*log(iSP)
#    return(sum(loglklhd))
# }

#=sa.lt.ku=================================================================================================
#simulated annealing for four parameters vitality model
#R function, version lt.514

sa.lt.ku<-function(r0,s0,k0,u0,x1,x2,Ni){
   
   #step 0 (initialization)
   n.p<-4   #number of parameters
   x0<-c(r0,s0,k0,u0)   #starting point

   v<-c(0.1,0.1,0.005,0.1)      # starting step vector v0, ...waiting for giving values
   T0<-10^3   #starting temperature T0
   eps<-10^(-7)  #a terminating criterion eps
   Ne<-4   #values of minima are less than a tolerance
   Ns<-20  #a test for step variation 
   c<-2    #a varying criterion 
   Nt<-10*n.p  #a test for temperature reduction
   rt<-0.85   #a reduction coefficient rt

   f0<-logLikelihood.ku(par=c(x0[1],x0[2],x0[3],x0[4]),x1,x2,Ni)  ###check the function form in splus
   x.opt<-x0
   x.temp<-x0
   x.new<-x0
   f.opt<-f0
   f.new<-f0
   f.temp<-f0
   N<-rep(0,n.p)
   x.star<-matrix(0,10000,n.p)
   x.star[1:Ne,]<-t(matrix(rep(x0,Ne),4,Ne))
   f.star<-c(rep(f0,Ne),rep(0,99996))

   i<-0   #successive points 
   j<-0   #successive cycles along every direction
   m<-0   #successive step adjustments
   k<-0   #successive temperature reductions
   h<-1   #the direction along which the trail point is generated

###############################################################
   #step 1
   while(k<1000){
            while(m<Nt){
                  while (j<Ns){
                        while (h<=n.p){
                              x.new[h]<-x.temp[h]+runif(1,-1,1)*v[h]
                              if (h==1){
                                   while(x.new[h]<0|x.new[h]>40){ #make sure all the points are in their suitable ranges
                                        x.new[h]<-x.temp[h]+runif(1,-1,1)*v[h] #generate a new point
                                   }
                              }
                              if (h==2){
                                   while(x.new[h]<0|x.new[h]>5){ #make sure all the points are in their suitable ranges
                                        x.new[h]<-x.temp[h]+runif(1,-1,1)*v[h] #generate a new point
                                   }
                              }
                              if (h==3){
                                   while(x.new[h]<0|x.new[h]>5){ #make sure all the points are in their suitable ranges
                                         x.new[h]<-x.temp[h]+runif(1,-1,1)*v[h] #generate a new point
                                   }
                              }
                              if (h==4){
                                   while(x.new[h]<0|x.new[h]>2){ #make sure all the points are in their suitable ranges
                                         x.new[h]<-x.temp[h]+runif(1,-1,1)*v[h] #generate a new point
                                   }
                              }
                             f.new<-logLikelihood.ku(par=c(x.new[1],x.new[2],x.new[3],x.new[4]),x1,x2,Ni)
     
                             if(f.new<=f.temp){         #then accept the new point
                                  x.temp<-x.new
                                  f.temp<-f.new
                                  i<-i+1
                                  N[h]<-N[h]+1
                                  if(f.new<f.opt){
                                      x.opt<-x.new
                                      f.opt<-f.new
                                  }
                             }
                            else{                      #metropolis move
                                 p<-exp((f.temp-f.new)/T0) 
                                 if(runif(1)<p){         #accept point
                                     x.temp<-x.new
                                     f.temp<-f.new
                                     i<-i+1
                                     N[h]<-N[h]+1
                                 }
                            }
                            h<-h+1
                        }
                        h<-1
                        j<-j+1
                 }
                 #for (t in 1:4){
                 #     if(N[t]>0.6*Ns){
                 #        v[t]<-v[t]*(1+c*(N[t]/Ns-0.6)/0.4)
                 #     }
                 #     else{
                 #         if(N[t]<0.4*NS){
                 #            v[t]<-v[t]/(1+c*(0.4-N[t]/Ns)/0.4)
                 #         }
                 #     }
                 #}
                j<-0
                N<-rep(0,n.p)
                m<-m+1
           }
           T0<-T0*rt
           f.star[k+Ne+1]<-f.temp
           x.star[k+Ne+1,]<-x.temp
           index<-0
           for (q in 1:Ne){
                if(abs(f.star[k+Ne+1]-f.star[k+Ne+1-q])>eps)
                   index<-1
           }
           for(t in 1:n.p){
                if(abs(x.star[k+Ne+1,t]-x.star[k+Ne,t])/x.star[k+Ne,t]>0.001)
                   index<-1
           }  
           if(f.star[k+Ne+1]-f.opt>eps)   
               index<-1
           if(index==0)  
               break                #stop the search 
           else{
               k<-k+1
               m<-0
           }
    }
    return(c(x.opt,f.opt))
}

#=stdErr.ku====================================================================================
#computing standand errors.

stdErr.ku<-function(r,s,k,u,x1,x2,Ni,pop){	
	# function to compute standard error for MLE parameters r,s,k,u in the vitality model
	# Arguments:
	#  r,s,k,u - final values of MLE parameters
	#  x1,x2  time vectors (steps 1:(T-1) and 2:T)
	#  Ni - survival fraction
	#  pop - total population of the study
	#
	# Return:
	#  standard error for r,s,k,u
	# Note: if k <or= 0, can not find std Err for k.
	#
	LL <-function(a,b,c,d,r,s,k,u,x1,x2,Ni){logLikelihood.ku(c(r+a,s+b,k+c,u+d),x1,x2,Ni)}
	
	#initialize hessian for storage
    hess <-matrix(0,nrow=4,ncol=4)

	
	#set finite difference intervals
	h <-.001
	hr <-abs(h*r)
	hs <-h*s*.1
	hk <-h*k*.1
	hu <-h*u*.1


	
	#Compute second derivitives (using 5 point)
	# LLrr
	f0 <-LL(-2*hr,0,0,0,r,s,k,u,x1,x2,Ni)
	f1 <-LL(-hr,0,0,0,r,s,k,u,x1,x2,Ni)
	f2 <-LL(0,0,0,0,r,s,k,u,x1,x2,Ni)
	f3 <-LL(hr,0,0,0,r,s,k,u,x1,x2,Ni)
	f4 <-LL(2*hr,0,0,0,r,s,k,u,x1,x2,Ni)

	fp0 <-(-25*f0 +48*f1 -36*f2 +16*f3 -3*f4)/(12*hr)
	fp1 <-(-3*f0 -10*f1 +18*f2 -6*f3 +f4)/(12*hr)
	fp3 <-(-f0 +6*f1 -18*f2 +10*f3 +3*f4)/(12*hr)
	fp4 <-(3*f0 -16*f1 +36*f2 -48*f3 +25*f4)/(12*hr)
	
	LLrr <-(fp0 -8*fp1 +8*fp3 -fp4)/(12*hr)

	# LLss
	f0 <-LL(0,-2*hs,0,0,r,s,k,u,x1,x2,Ni)
	f1 <-LL(0,-hs,0,0,r,s,k,u,x1,x2,Ni)
	# f2 as above
	f3 <-LL(0,hs,0,0,r,s,k,u,x1,x2,Ni)
	f4 <-LL(0,2*hs,0,0,r,s,k,u,x1,x2,Ni)

	fp0 <-(-25*f0 +48*f1 -36*f2 +16*f3 -3*f4)/(12*hs)
	fp1 <-(-3*f0 -10*f1 +18*f2 -6*f3 +f4)/(12*hs)
	fp3 <-(-f0 +6*f1 -18*f2 +10*f3 +3*f4)/(12*hs)
	fp4 <-(3*f0 -16*f1 +36*f2 -48*f3 +25*f4)/(12*hs)
	
	LLss <-(fp0 -8*fp1 +8*fp3 -fp4)/(12*hs)
	
	# LLkk

		f0 <-LL(0,0,-2*hk,0,r,s,k,u,x1,x2,Ni)
		f1 <-LL(0,0,-hk,0,r,s,k,u,x1,x2,Ni)
		# f2 as above
		f3 <-LL(0,0,hk,0,r,s,k,u,x1,x2,Ni)
		f4 <-LL(0,0,2*hk,0,r,s,k,u,x1,x2,Ni)

		fp0 <-(-25*f0 +48*f1 -36*f2 +16*f3 -3*f4)/(12*hk)
		fp1 <-(-3*f0 -10*f1 +18*f2 -6*f3 +f4)/(12*hk)
		fp3 <-(-f0 +6*f1 -18*f2 +10*f3 +3*f4)/(12*hk)
		fp4 <-(3*f0 -16*f1 +36*f2 -48*f3 +25*f4)/(12*hk)
	
		LLkk <-(fp0 -8*fp1 +8*fp3 -fp4)/(12*hk)
	
	
	# LLuu
	f0 <-LL(0,0,0,-2*hu,r,s,k,u,x1,x2,Ni)
	f1 <-LL(0,0,0,-hu,r,s,k,u,x1,x2,Ni)
	# f2 as above
	f3 <-LL(0,0,0,hu,r,s,k,u,x1,x2,Ni)
	f4 <-LL(0,0,0,2*hu,r,s,k,u,x1,x2,Ni)

	fp0 <-(-25*f0 +48*f1 -36*f2 +16*f3 -3*f4)/(12*hu)
	fp1 <-(-3*f0 -10*f1 +18*f2 -6*f3 +f4)/(12*hu)
	fp3 <-(-f0 +6*f1 -18*f2 +10*f3 +3*f4)/(12*hu)
	fp4 <-(3*f0 -16*f1 +36*f2 -48*f3 +25*f4)/(12*hu)
	
	LLuu <-(fp0 -8*fp1 +8*fp3 -fp4)/(12*hu)
	
	#-------end second derivs---
	# do mixed partials (4 points)
	# LLrs
	m1 <-LL(hr,hs,0,0,r,s,k,u,x1,x2,Ni)
	m2 <-LL(-hr,hs,0,0,r,s,k,u,x1,x2,Ni)
	m3 <-LL(-hr,-hs,0,0,r,s,k,u,x1,x2,Ni)
	m4 <-LL(hr,-hs,0,0,r,s,k,u,x1,x2,Ni)

	LLrs <-(m1 -m2 +m3 -m4)/(4*hr*hs)

   # LLru
	m1 <-LL(hr,0,0,hu,r,s,k,u,x1,x2,Ni)
	m2 <-LL(-hr,0,0,hu,r,s,k,u,x1,x2,Ni)
	m3 <-LL(-hr,0,0,-hu,r,s,k,u,x1,x2,Ni)
	m4 <-LL(hr,0,0,-hu,r,s,k,u,x1,x2,Ni)

	LLru <-(m1 -m2 +m3 -m4)/(4*hr*hu)
	
	# LLsu
	m1 <-LL(0,hs,0,hu,r,s,k,u,x1,x2,Ni)
	m2 <-LL(0,-hs,0,hu,r,s,k,u,x1,x2,Ni)
	m3 <-LL(0,-hs,0,-hu,r,s,k,u,x1,x2,Ni)
	m4 <-LL(0,hs,0,-hu,r,s,k,u,x1,x2,Ni)

	LLsu <-(m1 -m2 +m3 -m4)/(4*hu*hs)
	
	
		# LLrk
		m1 <-LL(hr,0,hk,0,r,s,k,u,x1,x2,Ni)
		m2 <-LL(-hr,0,hk,0,r,s,k,u,x1,x2,Ni)
		m3 <-LL(-hr,0,-hk,0,r,s,k,u,x1,x2,Ni)
		m4 <-LL(hr,0,-hk,0,r,s,k,u,x1,x2,Ni)

		LLrk <-(m1 -m2 +m3 -m4)/(4*hr*hk)
	
		# LLsk
		m1 <-LL(0,hs,hk,0,r,s,k,u,x1,x2,Ni)
		m2 <-LL(0,-hs,hk,0,r,s,k,u,x1,x2,Ni)
		m3 <-LL(0,-hs,-hk,0,r,s,k,u,x1,x2,Ni)
		m4 <-LL(0,hs,-hk,0,r,s,k,u,x1,x2,Ni)

		LLsk <-(m1 -m2 +m3 -m4)/(4*hs*hk)
		
		# LLku
		m1 <-LL(0,0,hk,hu,r,s,k,u,x1,x2,Ni)
		m2 <-LL(0,0,hk,-hu,r,s,k,u,x1,x2,Ni)
		m3 <-LL(0,0,-hk,-hu,r,s,k,u,x1,x2,Ni)
		m4 <-LL(0,0,-hk,hu,r,s,k,u,x1,x2,Ni)

		LLku <-(m1 -m2 +m3 -m4)/(4*hu*hk)


	
		diag(hess) <-c(LLrr,LLss,LLkk,LLuu)*pop
		hess[2,1]=hess[1,2]<-LLrs*pop
		hess[3,1]=hess[1,3]<-LLrk*pop
		hess[3,2]=hess[2,3]<-LLsk*pop
		hess[4,1]=hess[1,4]<-LLru*pop
		hess[4,2]=hess[2,4]<-LLsu*pop
		hess[4,3]=hess[3,4]<-LLku*pop


	#print(hess)
	hessInv <-solve(hess)
	#print(hessInv)

   #compute correlation matrix:
	sz <-4
	corr<-matrix(0,nrow=sz,ncol=sz)
	for (i in 1:sz) {
		for (j in 1:sz) {
			corr[i,j] <-hessInv[i,j]/sqrt(abs(hessInv[i,i]*hessInv[j,j]))
		}
	}
	#print(corr)
	if ( abs(corr[2,1]) > .98 ) {
		warning("WARNING: parameters r and s appear to be closely correlated for this data set.  
		   s.e. may fail for these parameters.")
	}
	if (  sz == 4 && abs(corr[3,2]) > .98 ) {
		warning("WARNING: parameters s and k appear to be closely correlated for this data set.  
		   s.e. may fail for these parameters.")
	}
	if (  sz == 4 && abs(corr[3,1]) > .98 ) {
		warning("WARNING: parameters r and k appear to be closely correlated for this data set.  
		   s.e. may fail for these parameters.")
	}
	if (  sz == 4 && abs(corr[4,2]) > .98 ) {
		warning("WARNING: parameters s and u appear to be closely correlated for this data set.  
		   s.e. may fail for these parameters.")
	}
	if (  sz == 4 && abs(corr[4,1]) > .98 ) {
		warning("WARNING: parameters r and u appear to be closely correlated for this data set.  
		   s.e. may fail for these parameters.")
	}
   
   	se <-sqrt(diag(hessInv))
	
    #  Approximate s.e. for cases where calculation of s.e. failed:
	 if( sum( is.na(se) ) > 0 ) {
             seNA<-is.na(se)
		 se12 <-sqrt(diag(solve(hess[c(1,2)	,c(1,2) ])))
		 se13 <-sqrt(diag(solve(hess[c(1,3)	,c(1,3) ])))
		 se23 <-sqrt(diag(solve(hess[c(2,3)	,c(2,3) ])))
		 se14 <-sqrt(diag(solve(hess[c(1,4)	,c(1,4) ])))
     se24 <-sqrt(diag(solve(hess[c(2,4)	,c(2,4) ])))
     se34 <-sqrt(diag(solve(hess[c(3,4)	,c(3,4) ])))
		
		 if(seNA[1]) {
			if(!is.na(se12[1]) ){
                     se[1]=se12[1]
                     warning("* s.e. for parameter r is approximate.")
                  }
                  else if(!is.na(se13[1])){
                          se[1]=se13[1]
                          warning("* s.e. for parameter r is approximate.")
                       }   
                       else if(!is.na(se14[1])){
                               se[1]=se14[1]
                               warning("* s.e. for parameter r is approximate.")
                            }   
				    else warning("* unable to calculate or approximate s.e. for parameter r.")
		 }
		 if(seNA[2]) {
			if(!is.na(se12[2]) ){
                     se[2]=se12[2]
                     warning("* s.e. for parameter s is approximate.")
                  }
                  else if(!is.na(se23[1])){
                          se[2]=se23[1]
                          warning("* s.e. for parameter s is approximate.")
                       }   
                       else if(!is.na(se24[1])){
                               se[2]=se24[1]
                               warning("* s.e. for parameter s is approximate.")
                            }   
				    else warning("* unable to calculate or approximate s.e. for parameter s.")
		 }
		 if(seNA[3]) {
			if(!is.na(se13[2]) ){
                     se[3]=se13[2]
                     warning("* s.e. for parameter k is approximate.")
                  }
                  else if(!is.na(se23[2])){
                          se[3]=se23[2]
                          warning("* s.e. for parameter k is approximate.")
                       }   
                       else if(!is.na(se34[1])){
                               se[3]=se34[1]
                               warning("* s.e. for parameter k is approximate.")
                            }   
				    else warning("* unable to calculate or approximate s.e. for parameter k.")
		 }
             if(seNA[4]) {
			if(!is.na(se14[2]) ){
                     se[4]=se14[2]
                     warning("* s.e. for parameter u is approximate.")
                  }
                  else if(!is.na(se24[2])){
                          se[4]=se24[2]
                          warning("* s.e. for parameter u is approximate.")
                       }   
                       else if(!is.na(se34[1])){
                               se[4]=se34[2]
                               warning("* s.e. for parameter u is approximate.")
                            }   
				    else warning("* unable to calculate or approximate s.e. for parameter u.")
		 }

	}


   	#######################
	return(se)
}

#=plotting.ku=============================================================================================
plotting.ku<-function(r.final,s.final,k.final,u.final,mlv,time,sfract,x1,x2,Ni,pplot,tlab,lplot,cplot,Iplot,gfit){
	#  Function to provide plotting and goodness of fit computations
	
# --plot cumulative survival---
	if(pplot != F){
    #win.graph()
		ext<-max(pplot,1)
   		par(mfrow=c(1,1))
   		len<-length(time)
   		tmax <-ext*time[len]
   		plot(time,sfract,xlab=tlab,ylab="survival fraction",ylim=c(0,1),xlim=c(0,tmax))
		xxx<-seq(0,tmax,length=200)
   		lines(xxx,SurvFn.ku(xxx,r.final,s.final,k.final,u.final))
   		title("Cumulative Survival Data and Vitality Model Fitting")

       # --likelihood and likelihood contour plots---
		if(lplot != F) {
		  profilePlot.ku<-function(r.f,s.f,k.f,u.f,x1,x2,Ni,mlv,cplot) {
		    #
		    # mlv = value of max likelihood 
		    # likelihood plots
		    # These plots are really not "profile" plots in that while one parameter is varied, 
		    # the others are held fixed, not re-optimized. 
		    
		    SLL <- function(r,s,k,u,x1,x2,Ni){logLikelihood.ku(c(r,s,k,u),x1,x2,Ni)}
		    
		    rf<-.2; sf<-.5; kf<-1.0; uf<-0.8;fp<-40 # rf,sf,kf,uf - set profile plot range (.2 => plot +-20%), 2*fp+1 points
		    rseq <-seq((1-rf)*r.f,(1+rf)*r.f, (rf/fp)*r.f)
		    sseq <-seq((1-sf)*s.f,(1+sf)*s.f, (sf/fp)*s.f)
		    useq <-seq((1-uf)*u.f,(1+uf)*u.f, (uf/fp)*u.f)
		    
		    if (k.f > 0) {
		      kseq <-seq((1-kf)*k.f,(1+kf)*k.f, (kf/fp)*k.f)
		    }
		    else {  #if k=0..
		      kseq <-seq(.00000001,.1,length=(2*fp+1))
		    }
		    
		    rl <-length(rseq)
		    tmpLLr <-rep(0,rl)
		    tmpLLs <-tmpLLr
		    tmpLLk <-tmpLLr
		    tmpLLu <-tmpLLr
		    for (i in 1:rl) {
		      tmpLLr[i] <-SLL(rseq[i],s.f,k.f,u.f,x1,x2,Ni)
		      tmpLLs[i] <-SLL(r.f,sseq[i],k.f,u.f,x1,x2,Ni)	
		      tmpLLk[i] <-SLL(r.f,s.f,kseq[i],u.f,x1,x2,Ni)
		      tmpLLu[i] <-SLL(r.f,s.f,k.f,useq[i],x1,x2,Ni)
		    }
		    
		    par(mfrow=c(2,2), mar=c(5,4,3,2))
		    
		    rlim1 <-rseq[1]
		    rlim2 <-rseq[rl]
		    if (r.f < 0) {  #even though r should not be <0
		      rlim2 <-rseq[1]
		      rlim1 <-rseq[rl]
		    }
		    LL<-SLL(r.f,s.f,k.f,u.f,x1,x2,Ni)		  
		    plot(r.f,LL, xlim=c(rlim1,rlim2), xlab="r",ylab="Likelihood")  
		    title("Likelihood Plots", outer=T, line=-1)
		    lines(rseq,tmpLLr)
		    legend(x="topright", legend=c("r.final", "vary r"), pch=c(1,NA), lty=c(NA,1))
		    plot(s.f,LL,xlim=c(sseq[1],sseq[rl]), xlab="s",ylab="Likelihood") 
		    lines(sseq,tmpLLs)
		    legend(x="topright", legend=c("s.final", "vary s"), pch=c(1,NA), lty=c(NA,1))
		    plot(k.f,LL,xlim=c(kseq[1],kseq[rl]), ylim=c(1.1*min(tmpLLk)-.1*(mLk=max(tmpLLk)),mLk),xlab="k",ylab="Likelihood") 
		    lines(kseq,tmpLLk)
		    legend(x="topright", legend=c("k.final", "vary k"), pch=c(1,NA), lty=c(NA,1))
		    plot(u.f,LL,xlim=c(useq[1],useq[rl]), xlab="u",ylab="Likelihood")
		    lines(useq,tmpLLu)
		    legend(x="topright", legend=c("u.final", "vary u"), pch=c(1,NA), lty=c(NA,1))
		    
		    
		    # --  for contour plotting  ----------------------------------------
		    
		    if (cplot==T) {
		      rl2<-(rl+1)/2; rl4<-20; st<-rl2-rl4; nr<-2*rl4+1
		      tmpLLrs <-matrix(rep(0,nr*nr),nrow<-nr,ncol=nr)
		      tmpLLru <-matrix(rep(0,nr*nr),nrow<-nr,ncol=nr)
		      tmpLLsu <-matrix(rep(0,nr*nr),nrow<-nr,ncol=nr)
		      for(i in 1:nr) {
		        for(j in 1:nr) {
		          tmpLLrs[i,j] <-SLL(rseq[i+st-1],sseq[j+st-1],k.f,u.f,x1,x2,Ni)
		        }
		      }
		      
		      for(i in 1:nr) {
		        for(j in 1:nr) {
		          tmpLLru[i,j] <-SLL(rseq[i+st-1],s.f,k.f,useq[j+st-1],x1,x2,Ni)
		        }
		      }	
		      for(i in 1:nr) {
		        for(j in 1:nr) {
		          tmpLLsu[i,j] <-SLL(r.f,sseq[i+st-1],k.f,useq[j+st-1],x1,x2,Ni)
		        }
		      }	
		      lvv<-seq(mlv,1.02*mlv,length=11)  #99.8%, 99.6% ... 98%
		      par(mfrow=c(2,2))
		      contour(rseq[st:(rl2+rl4)],sseq[st:(rl2+rl4)],tmpLLrs,levels=lvv,xlab="r",ylab="s")
		      title("Likelihood Contour Plot of r and s.")
		      points(r.f,s.f,pch="*",cex=3.0)
		      points(c(rseq[st],r.f),c(s.f,sseq[st]),pch="+",cex=1.5)
		      
		      contour(rseq[st:(rl2+rl4)],useq[st:(rl2+rl4)],tmpLLru,levels=lvv,xlab="r",ylab="u")
		      title("Likelihood Contour Plot of r and u.")
		      points(r.f,u.f,pch="*",cex=3.0)
		      points(c(rseq[st],r.f),c(u.f,useq[st]),pch="+",cex=1.5)
		      
		      contour(sseq[st:(rl2+rl4)],useq[st:(rl2+rl4)],tmpLLsu,levels=lvv,xlab="s",ylab="u")
		      title("Likelihood Contour Plot of s and u.")
		      points(s.f,u.f,pch="*",cex=3.0)
		      points(c(sseq[st],s.f),c(u.f,useq[st]),pch="+",cex=1.5)
		      
		      plot(1,1, pch=NA,xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
		      legend(x="center", legend="Outermost ring is likelihood 98% (of max) level,\ninnermost is 99.8% level.", bty="n")
		    }
		    
		  }
			profilePlot.ku(r.final,s.final,k.final,u.final,x1,x2,Ni,mlv,cplot)
		}
	}

# --calculations for goodness of fit---
	if(gfit!=F){ # then gfit must supply the population number
   		isp <-survProbInc.ku(r.final,s.final,k.final,u.final,x1,x2)
   		C1.calc<-function(pop, isp, Ni){
   		  #  Routine to calculate goodness of fit (Pearson's C -type test)
   		  #  pop - population number of sample
   		  #  isp  - Modeled: incemental survivor Probability 
   		  #  Ni   - Data:  incemental survivor fraction (prob)
   		  #
   		  # Returns:  a list containing:
   		  # C1, dof, Chi2    (retrieve each from list as ..$C1   etc.)
   		  
   		  if (pop < 35) {
   		    if (pop < 25) {
   		      warning(paste("WARNING: sample population (",as.character(pop),") is too small for 
   		                    meaningful goodness of fit measure.  Goodness of fit not being computed"))
   		      return()
   		    } else {
   		      warning(paste("WARNING: sample population (",as.character(pop),") may be too small for 
   		                    meaningful goodness of fit measure"))
   		    }
   		    }
   		  np <- pop * isp # modeled population at each survival probability level.
   		  tmpC1 <- 0
   		  i1 <-1; i<-1; cnt<-0
   		  len <- length(np)
   		  while(i <= len) {
   		    idx <- i1:i
   		    # It is recommended that each np[i] >5 for meningful results.  Where np<5
   		    #   points are grouped to attain that level.   I have fudged it to 4.5 ... 
   		    #  (as some leeway is allowed, and exact populations are sometimes unknown).
   		    if(sum(np[idx]) > 4.5) {
   		      cnt <-cnt+1
   		      # Check if enough points remain.  If not, they are glommed onto previous grouping.
   		      if(i < len && sum(np[(i + 1):len]) < 4.5) {
   		        idx <- i1:len
   		        i <- len
   		      }
   		      
   		      sNi <- sum(Ni[idx])
   		      sisp <- sum(isp[idx])
   		      tmpC1 <- tmpC1 + (pop * (sNi - sisp)^2)/sisp
   		      i1 <-i+1			
   		    }
   		    
   		    i <-i+1
   		  }
   		  C1 <- tmpC1
   		  dof <-cnt-1-4    # degrees of freedom (3 is number of parameters).
   		  if (dof < 1) {
   		    warning(paste("WARNING: sample population (",as.character(pop),") is too small for 
		meaningful goodness of fit measure (DoF<1).  Goodness of fit not being computed"))
   		    return()
   		  }
   		  
   		  chi2<-qchisq(.95,dof)
   		  return(list(C1=C1,dof=dof,chi2=chi2))
   		  }
		C1dof <- C1.calc(gfit,isp,Ni)
		C1 <-C1dof$C1
		dof <-C1dof$dof
		chi2 <-C1dof$chi2
				
		print(paste("Pearson's C1=",as.character(round(C1,3)),"   chisquared =",
				as.character(round(chi2,3)),"on",as.character(dof),"degrees of freedom "))
		
		#   Note: The hypothesis being tested is whether the data could reasonably have come 
		#   from the assumed (vitality) model.
		if (C1 > chi2) {
			print("C1 > chiSquared;   should reject the hypothesis becasue C1 falls outside the 95% confidence interval.")
		} else {
			print("C1 < chiSquared;   should Not reject the hypothesis because C1 falls inside the 95% confidence interval")		
		}
	}
		
	# --Incremental mortality plot
	if(Iplot != F){
		par(mfrow=c(1,1))
		#if (rc.data != F) {
			ln <-length(Ni)-1
			x1 <-x1[1:ln]
			x2 <-x2[1:ln]
			Ni <-Ni[1:ln]
		#}
		#ln <-length(Ni)
		#scale<-(x2-x1)[Ni==max(Ni)]
		scale<-max( (x2-x1)[Ni==max(Ni)] )
					
		ext<-max(pplot,1)
		
		npt<-200*ext
		xxx <-seq(x1[1],x2[ln]*ext,length=npt)
		xx1 <-xxx[1:(npt-1)]
		xx2 <-xxx[2:npt]
		sProbI <-survProbInc.ku(r.final,s.final,k.final,u.final,xx1,xx2)
			
		ytop <-1.1*max( max(sProbI/(xx2-xx1)),Ni/(x2-x1) )*scale		
		plot((x1+x2)/2,Ni*scale/(x2-x1),ylim=c(0,ytop),xlim=c(0,ext*x2[ln]),xlab=tlab,ylab="incremental mortality")
		title("Probability Density Function")
		lines((xx1+xx2)/2,sProbI*scale/(xx2-xx1))
	}
	
	return()	 
}
	
