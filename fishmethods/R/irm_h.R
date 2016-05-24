irm_h<-function(relyrs=NULL,recapyrs=NULL,N=NULL,recapharv=NULL,lambda=NULL,phi=NULL, 
        Fyr=NULL, Myr=NULL, initial=NULL,lower=c(0.0001,0.0001),upper=c(5,5),maxiter=10000){
      Fp<-rep(initial[1],length(Fyr))
      Mp<-rep(initial[2],length(Myr))  
### Error Messages
    if(class(recapharv)!="matrix") stop("recapharv is not a matrix.")
    if(is.null(relyrs)|is.null(recapyrs)) stop("Missing relyrs or recapyrs.")
    if(is.null(N)) stop("Ns are missing.")
    if(is.null(recapharv)) stop("Missing recovery matrix for harvested fish.")
    if(is.null(lambda)) stop("lambdas (reporting rates) for harvested fish) are missing.")
    if(is.null(phi)) stop("hphi (initial tag survival rates) for harvested fish) are missing.")
    if(is.null(Fyr)) stop("Year designations for fishing mortality estimates (Fyr) are missing.")
    if(is.null(Myr)) stop("Year designations for natural mortality estimates (Myr) are missing.")
    if(Myr[1]!=relyrs[1]|Myr[1]!=recapyrs[1]) stop("First year in Myr should be equal to the first year in relyrs and recapyrs.")
    if(Fyr[1]!=relyrs[1]|Fyr[1]!=recapyrs[1]) stop("First year in Fyr should be equal to the first year in relyrs and recapyrs.")
   if(relyrs[1]!=recapyrs[1]) stop("First year in relyrs and recapyrs should be the same.")
    rec<-length(seq(recapyrs[1],recapyrs[2],1))
    if(ncol(recapharv)!=rec) stop("The number of columns in recapharv does not equal the number of
             year specified by recapyrs.")
    if(length(lambda)!=rec) stop("The number of values in lambda does not equal the number of
             year specified by recapyrs.")
    if(length(phi)!=rec) stop("The number of values in phi does not equal the number of
             year specified by recapyrs.")
     rel<-length(seq(relyrs[1],relyrs[2],1))
      if(length(N)!=rel) stop("The number of values in N does not equal the number of
             year specified by relyrs.")
	if(initial[1]<lower[1]|initial[1]>upper[1]) stop("initial F is outside lower or upper bounds.")
	if(initial[2]<lower[2]|initial[2]>upper[2]) stop("initial M is outside lower or upper bounds.")
      
##########Program
      parms<-c(Fp,Mp)
     #x<-parms
   fit.obs<-function(x){
        F<-x[1:length(Fp)]
        M<-x[as.numeric(length(Fp)+1):as.numeric(length(Fp)+length(Mp))]
#######calc_number_tags
      	styrR<-1;endyrR<-relyrs[2]-relyrs[1]+1
      	styr<-1;endyr<-recapyrs[2]-recapyrs[1]+1
     		tags<-NULL
     		cnt<-0
    		 for (t in styrR:endyrR){
    	 	 Ntags<-0
   		 for (y in as.numeric(styr+cnt):endyr){
       	     Ntags<-Ntags+recapharv[t,y]      	   }
    		   tags[t]<-Ntags
       	   cnt<-cnt+1
       	}
     		yrvector<-seq(recapyrs[1],recapyrs[2],1)    
#####calc_F_vector
  		Fvector<-rep(NA,length(yrvector)) 
  		Ftemp<-Fyr
  		Ftemp[length(Ftemp)+1]<-recapyrs[2]+1
 		for(t in styr:endyr){
     			for(d in 1:length(F)){
         		if(yrvector[t]>=Ftemp[d] & yrvector[t]<Ftemp[d+1]) Fvector[t]<-F[d]
       		}
   		}
###### calc_M_vector
  		Mvector<-rep(NA,length(yrvector)) 
  		Mtemp<-Myr
  		Mtemp[length(Mtemp)+1]<-recapyrs[2]+1
 		for(t in styr:endyr){
     			for(d in 1:length(M)){
         		if(yrvector[t]>=Mtemp[d] & yrvector[t]<Mtemp[d+1]) Mvector[t]<-M[d]
       		}
   		}
####calc_s
  		s<-array(NA,dim=c(endyrR,endyr))
   		cnt<-0
  		for (t in styrR:endyrR){
    			for (y in as.numeric(styr+cnt):endyr){
       			if(t==y){s[t,y]<-1}
       			if(t!=y){
         			s[t,y]<-exp(-Fvector[y-1]-Mvector[y-1])    
        			}
     		 	}		   
      		cnt<-cnt+1
   		}
####calc_u_h
  		u_h<-array(NA,dim=c(endyrR,endyr))
  		cnt<-0
 		for (t in styrR:endyrR){
   			for (y in as.numeric(styr+cnt):endyr){
     	  			u_h[t,y]<-(Fvector[y]/(Fvector[y]+Mvector[y]))*(1-exp(-Fvector[y]-Mvector[y]))
      		}		   
      		cnt<-cnt+1
   		}
####calc_s_prob
 		cnt<-0
 		s_prob<-array(NA,dim=c(endyrR,endyr))
 		for (t in styrR:endyrR){
    		looper<-0
    			for (y in as.numeric(styr+cnt):endyr){
				probs<-1
				for(a in as.numeric(y-looper):y){
           				probs<-probs*s[t,a]
          			}
          		s_prob[t,y]<-probs
          		looper<-looper+1
      		}		   
    			cnt<-cnt+1
   		}
####calc_exp_prob_h
 		exp_prob_h<-array(NA,dim=c(endyrR,endyr))
 		sum_prob_h<-NULL
  		cnt<-0
 		for (t in styrR:endyrR){
    			dodo<-0
    			for (y in as.numeric(styr+cnt):endyr){
        			exp_prob_h[t,y]<-lambda[y]*phi[y]*s_prob[t,y]*u_h[t,y]
	  			dodo<-dodo+exp_prob_h[t,y]
      		}	
    			sum_prob_h[t]<-dodo   
    			cnt<-cnt+1
   		} 
####calc_LL Likelihood calculation
  		ll_h<-array(NA,dim=c(endyrR,endyr))
  		ll_ns<-NULL
 		cnt<-0
 		for (t in styrR:endyrR){
    			for (y in as.numeric(styr+cnt):endyr){
       			ll_h[t,y]<-0
        			if(recapharv[t,y]!=0){
          				ll_h[t,y]<-recapharv[t,y]*log(exp_prob_h[t,y])
         			}
        		
      		}		   
    			cnt<-cnt+1
   		}
 		for (t in styrR:endyrR){
     			ll_ns[t]<-(N[t]-tags[t])*log(1-sum_prob_h[t])
   		}
####Sum likelihood
 		f_tag<-0
 		cnt<-0
 		for (t in styrR:endyrR){
    			for (y in as.numeric(styr+cnt):endyr){
       			f_tag<-f_tag+ll_h[t,y]  
      		}		 
    			cnt<-cnt+1
   		}
 		for (t in styrR:endyrR){
       		f_tag<-f_tag+ll_ns[t]
   		} 
    		LL<-f_tag*-1
    		LL
 	}#end fit.obs function

########################Estimate
	low<-c(rep(lower[1],length(Fyr)),rep(lower[2],length(Myr)))
	up<-c(rep(upper[1],length(Fyr)),rep(upper[2],length(Myr)))
	results<-optim(parms, fit.obs, gr = NULL,lower=low,upper=up,method=c("L-BFGS-B"), 
      control=list(maxit=maxiter),hessian=TRUE)
 	var<-diag(solve(results$hessian))
  	varcov<-solve(results$hessian)
 	cormat<-(cov2cor(varcov))
      corr<-list(c(paste("F",1:length(Fp),sep=""),
                   paste("M",1:length(Mp),sep="")))
     dimnames(cormat)[1]<-corr
     dimnames(cormat)[2]<-corr

############################CALCULATE OUTPUT##############################
   calc.outpt<-function(){

#######calc_number_tags
      styrR<-1;endyrR<-relyrs[2]-relyrs[1]+1
      styr<-1;endyr<-recapyrs[2]-recapyrs[1]+1
      tags<-NULL
      cnt<-0
      for (t in styrR:endyrR){
    	 Ntags<-0
   	 for (y in as.numeric(styr+cnt):endyr){
            Ntags<-Ntags+recapharv[t,y]
         }
       tags[t]<-Ntags
       cnt<-cnt+1
       }
######################Get F, FA, and M estimates and spread throughout years
 	F<-results$par[1:length(Fp)];FVAR<-var[1:length(Fp)]
 	M<-results$par[as.numeric(length(Fp)+1):as.numeric(length(Fp)+length(Mp))]
 	MVAR<-var[as.numeric(length(Fp)+1):as.numeric(length(Fp)+length(Mp))]
 	yrvector<-seq(recapyrs[1],recapyrs[2],1)  
      Fpos<-1:length(Fp)
      Mpos<-(length(Fp)+1):(length(Fp)+length(Mp))  
 #####calc_F_vector
  	Fvector<-data.frame(F=rep(NA,length(yrvector)),VAR=rep(NA,length(yrvector)),Fpos=rep(NA,length(yrvector))) 
  	Ftemp<-Fyr
  	Ftemp[length(Ftemp)+1]<-recapyrs[2]+1
 	for(t in styr:endyr){
     		for(d in 1:length(F)){
         		if(yrvector[t]>=Ftemp[d] & yrvector[t]<Ftemp[d+1]) {
              	Fvector[t,1]<-F[d]
              	Fvector[t,2]<-FVAR[d]
                  Fvector[t,3]<-Fpos[d]
          		}
       	}
   	}
###### calc_M_vector
  	Mvector<-data.frame(M=rep(NA,length(yrvector)),VAR=rep(NA,length(yrvector)),Mpos=rep(NA,length(yrvector))) 
  	Mtemp<-Myr
  	Mtemp[length(Mtemp)+1]<-recapyrs[2]+1
 	for(t in styr:endyr){
     		for(d in 1:length(M)){
         		if(yrvector[t]>=Mtemp[d] & yrvector[t]<Mtemp[d+1]){
            	Mvector[t,1]<-M[d]
            	Mvector[t,2]<-MVAR[d]
			Mvector[t,3]<-Mpos[d]
          		}
       	}
   	}
###########calculate S vector
 	Zvector<-cbind(Fvector,Mvector)
       names(Zvector)<-c("F","FVAR","Fpos","M","MVAR","Mpos")
       
 	for(i in 1:nrow(Zvector)){
    		Zvector$FM[i]<-varcov[Zvector$Fpos[i],Zvector$Mpos[i]]
  	}

  	for(i in 1:nrow(Zvector)){
 		Zvector$Z[i]<-Zvector$F[i]+Zvector$M[i]
 		Zvector$ZVAR[i]<-Zvector$FVAR[i]+Zvector$MVAR[i]+
                  2*Zvector$FM[i]
 	}
 	Zvector$S<-exp(-Zvector$Z)
 	Zvector$SVAR<-Zvector$ZVAR*exp(-Zvector$Z)^2
 	Svector<-data.frame(S=Zvector$S,VAR=Zvector$SVAR,SE=sqrt(Zvector$SVAR))
       Mvector$SE<-sqrt(Mvector$VAR)
       Mvector<-Mvector[,-c(3)]
       Mvector$Year<-seq(recapyrs[1],recapyrs[2],1)   
       Fvector$SE<-sqrt(Fvector$VAR)
       Fvector<-Fvector[,-c(3)]
       Fvector$Year<-seq(recapyrs[1],recapyrs[2],1)
       Zvector<-Zvector[,8:9]
       Zvector$SE<-sqrt(Zvector$ZVAR)
       names(Zvector)<-c("Z","VAR","SE")
	 Zvector$Year<-seq(recapyrs[1],recapyrs[2],1)
       Svector$Year<-seq(recapyrs[1],recapyrs[2],1)      
 ####calc_s
  	s<-array(NA,dim=c(endyrR,endyr))
   	cnt<-0
  	for (t in styrR:endyrR){
    		for (y in as.numeric(styr+cnt):endyr){
       		if(t==y){s[t,y]<-1}
       		if(t!=y){
         		s[t,y]<-exp(-Fvector[y-1,1]-Mvector[y-1,1])    
        		}
      	}		   
      	cnt<-cnt+1
   	}
####calc_u_h
  	u_h<-array(NA,dim=c(endyrR,endyr))
  	cnt<-0
 	for (t in styrR:endyrR){
   		for (y in as.numeric(styr+cnt):endyr){
       	u_h[t,y]<-(Fvector[y,1]/(Fvector[y,1]+Mvector[y,1]))*(1-exp(-Fvector[y,1]-Mvector[y,1]))
      	}		   
      	cnt<-cnt+1
   	}
####calc_s_prob
 	cnt<-0
 	s_prob<-array(NA,dim=c(endyrR,endyr))
 	for (t in styrR:endyrR){
    		looper<-0
    		for (y in as.numeric(styr+cnt):endyr){
			probs<-1
			for(a in as.numeric(y-looper):y){
           			probs<-probs*s[t,a]
          		}
          		s_prob[t,y]<-probs
          		looper<-looper+1
      	}		   
    		cnt<-cnt+1
   	}
####calc_exp_prob_h
 	exp_prob_h<-array(NA,dim=c(endyrR,endyr))
 	sum_prob_h<-NULL
  	cnt<-0
 	for (t in styrR:endyrR){
    		dodo<-0
    		for (y in as.numeric(styr+cnt):endyr){
        		exp_prob_h[t,y]<-lambda[y]*phi[y]*s_prob[t,y]*u_h[t,y]
	  		dodo<-dodo+exp_prob_h[t,y]
      	}	
    		sum_prob_h[t]<-dodo   
    		cnt<-cnt+1
   	} 
############DIAGNOSTICS#####################################
####calc_Chisquare
   	styrR<-1;endyrR<-relyrs[2]-relyrs[1]+1
      styr<-1;endyr<-recapyrs[2]-recapyrs[1]+1
	exp_r_h<-array(NA,dim=c(endyrR,endyr))
 	cnt<-0;up_count<-0
 	for (t in styrR:endyrR){
    		for (y in as.numeric(styr+cnt):endyr) {
       		up_count<-up_count+1 
      	}	
      	cnt<-cnt+1
   	}
 	cnt<-0
 	for (t in styrR:endyrR){
    		for (y in as.numeric(styr+cnt):endyr){
       		exp_r_h[t,y]<-exp_prob_h[t,y]*N[t]
      	}	
    		cnt<-cnt+1
   	}
  	cnt<-0
    	chi_h<-array(NA,dim=c(endyrR,endyr))
  	chi_ns<-NULL
  	pear_h<-array(NA,dim=c(endyrR,endyr))
  	pear_ns<-NULL
  	exp_ns<-NULL
  	for (t in styrR:endyrR){
    		for (y in as.numeric(styr+cnt):endyr){
        		chi_h[t,y]<-(recapharv[t,y]-exp_r_h[t,y])^2/exp_r_h[t,y]
        		pear_h[t,y]<-(recapharv[t,y]-exp_r_h[t,y])/sqrt(exp_r_h[t,y])
      	}	
     		cnt<-cnt+1
   	}
  	for (t in styrR:endyrR){
      	exp_ns[t]<-N[t]*(1-sum_prob_h[t])
   	}
  #Not seen chi
  	for (t in styrR:endyrR){
        chi_ns[t]<-0
        chi_ns[t]<-((N[t]-tags[t])-exp_ns[t])^2/exp_ns[t]
        pear_ns[t]<-((N[t]-tags[t])-exp_ns[t])/sqrt(exp_ns[t])
   	}
####total chi square
 	up_chi<-sum(chi_h,na.rm=T)+sum(chi_ns,na.rm=T)
 	K<-length(Fyr)+length(Myr)
 	up_df<-up_count*2-K;
 	up_chat<-up_chi/up_df
 	AIC<-2*results$value+2*K
 	AICc<-AIC+(2*K*(K+1))/(sum(N)-K-1)
#####calc_pooled_cells
# Pool harvested cells
   	pool_h_e<-array(NA,dim=c(endyrR,endyr))
   	pool_h<-array(NA,dim=c(endyrR,endyr))
  	cnt<-0
  	for (t in styrR:endyrR){ 
      	for(y in as.numeric(styr+cnt):endyr){
            	pool_h_e[t,y]<-0
            	pool_h[t,y]<-0
            	pool_h_e[t,y]<-exp_r_h[t,y]
            	pool_h[t,y]<-recapharv[t,y]   
        	}
      	cnt<-cnt+1
   	}
  	cnt<-0
  	hless<-0
  	for(t in styrR:endyrR){    
    		for(y in endyr:as.numeric(styr+cnt)){
          		if(pool_h_e[t,y]>=1){
             		pool_h[t,y]<-pool_h[t,y]
             		pool_h_e[t,y]<-pool_h_e[t,y]
            	}
          		if(pool_h_e[t,y]>=0 & pool_h_e[t,y]<1){ 
 		  		if (y!=as.numeric(styr+cnt)){
               		 	hless<-hless+1
                			pool_h_e[t,y-1]<-pool_h_e[t,y-1]+pool_h_e[t,y]
                			pool_h[t,y-1]<-pool_h[t,y-1]+pool_h[t,y]
                			pool_h[t,y]<-0
                			pool_h_e[t,y]<-0
               		}
               	if (y==as.numeric(styr+cnt)) next #change from break
            	}
         	}
         	cnt<-cnt+1
     }
   	p_df<-up_df
####Pooled Chi-square
   	p_chi_h<-array(NA,dim=c(endyrR,endyr))
   	cnt<-0
   	for (t in styrR:endyrR) {
    		for (y in as.numeric(styr+cnt):endyr){
        		p_chi_h[t,y]<-0
        		if(pool_h_e[t,y]!=0){
          			p_chi_h[t,y]<-(pool_h[t,y]-pool_h_e[t,y])^2/pool_h_e[t,y]
         		}
	     	}	
      	cnt<-cnt+1
	}
  	p_chi<-sum(p_chi_h,na.rm=T)+sum(chi_ns,na.rm=T)
  	p_chat<-p_chi/p_df
###########################Output##################################
     if(results$convergence==0) mess<-"Successful."
     if(results$convergence==1) mess<-"Maximum Iterations Reached."
     if(results$convergence==51) mess<-results$message
     if(results$convergence==53) mess<-results$message
	ans<-NULL
	ans$statistics<-matrix(NA,11L,1L)
	ans$statistics<-rbind(round(results$value*-1,3),round(K,0),round(AIC,2),round(AICc,2)
                  ,round(sum(N,0)),round(up_chi,2),round(up_df,0),round(up_chat,3),round(p_chi,2)
                  ,round(p_df,0),round(p_chat,3))   
	dimnames(ans$statistics)<-list(c("Neg. Log-Likelihood","K","AIC","AICc","Eff. Sample Size","Unpooled Chi-square",
              "Unpooled df","Unpooled c-hat","Pooled Chi-square","Pooled df","Pooled c-hat"))
	dimnames(ans$statistics)[[2]][1]<-list(c("Value"))
      ans$model_convergence<-mess
      ans$parameter_correlation_matrix<-cormat
	ans$fishing_mortality<-Fvector
      ans$natural_mortality<-Mvector
      ans$total_mortality<-Zvector
      ans$survival<-Svector
      colnames(recapharv)<-c(seq(recapyrs[1],recapyrs[2],1))
      rownames(recapharv)<-c(seq(relyrs[1],relyrs[2],1))
	ans$obs_recoveries_harvested<-recapharv
      colnames(exp_r_h)<-c(seq(recapyrs[1],recapyrs[2],1))
      rownames(exp_r_h)<-c(seq(relyrs[1],relyrs[2],1))
	ans$pred_recoveries_harvested<-exp_r_h
	ans$pred_number_notseen<-exp_ns	
      colnames(chi_h)<-c(seq(recapyrs[1],recapyrs[2],1))
      rownames(chi_h)<-c(seq(relyrs[1],relyrs[2],1))
	ans$unpooled_cell_chisquare_harvested<-chi_h
	ans$unpooled_cell_chisquare_notseen<-chi_ns
      colnames(pear_h)<-c(seq(recapyrs[1],recapyrs[2],1))
      rownames(pear_h)<-c(seq(relyrs[1],relyrs[2],1))
	ans$unpooled_cell_Pearson_harvested<-pear_h
	ans$unpooled_cell_Pearson_notseen<-pear_ns
      ans$type<-"h"
	return(ans)
  }#calc.outpt
 	calc.outpt()
}

