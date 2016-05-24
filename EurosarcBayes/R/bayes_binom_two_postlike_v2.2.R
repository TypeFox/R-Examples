# LINES trial exact "simulation" of expected outcome under specific data conditions.
# Since the number of patients is small a probabilistic approach is more efficient than a large sampling approach.
# Using probability of success to stop the trial early
# toxicity and efficacy considered independently
# note trial cannot stop if uncertainty over toxicity remains. This means that if there is high eficacy and toxicity then the design is likely to recommend stopping early due to high toxicity i.e. there is no built in trade-off within the model.



bayes_binom_two_postlike=function(t,r,reviews,pra=1,prb=1,pta=1,ptb=1,efficacy_critical_value=0.2,efficacy_prob_stop=0.9,toxicity_critical_value=0.3,toxicity_prob_stop=0.9,int_combined_prob=0,int_futility_prob=1,int_toxicity_prob=1,int_efficacy_prob=1,futility_critical_value=0.35,no_toxicity_critical_value=0.1){

  ##################################################
  # critical values at the end of the trial.
  n=max(reviews)

  R_cut=which(pbeta(efficacy_critical_value,pra+0:n,prb+n-0:n)<1-efficacy_prob_stop)[1]-1
  T_cut=max(which(pbeta(no_toxicity_critical_value,pta+0:n,ptb+n-0:n)>toxicity_prob_stop))-1

  # variable to store the critical values for each endpoint
  cutpoints=matrix(NA,length(reviews),5)
  cutpoints[,1]=reviews
  prob_r_sucess=list()
  prob_t_sucess=list()
  #################################################
  # interation free cutpoints and probability of sucess/fail under each category
  #################################################
  if(length(reviews)>1){
    for(i in 1:(length(reviews)-1)){
      prob_r_sucess[[i]]=pbetabinom.ab(R_cut-0:reviews[i]-1,n-reviews[i],pra+0:reviews[i],prb+reviews[i]-0:reviews[i])
      prob_r_sucess[[i]]=1-prob_r_sucess[[i]]
      for(j in 1:(reviews[i]+1)){
        if(prob_r_sucess[[i]][j]<0){
          prob_r_sucess[[i]][j]=0
        }
      }
      prob_t_sucess[[i]]=pbetabinom.ab(T_cut-0:reviews[i]-1,n-reviews[i],pta+0:reviews[i],ptb+reviews[i]-0:reviews[i])
      for(j in 1:(reviews[i]+1)){
        if(prob_t_sucess[[i]][j]<0){prob_t_sucess[[i]][j]=0}
      }

    }
  }



  decision=list()
  if(length(reviews)>1){
    for(i in 1:(length(reviews)-1)){
      decision[[i]]=matrix(0,reviews[i]+1,reviews[i]+1)
      for(xr in 1:(reviews[i]+1)){
        for(xt in 1:(reviews[i]+1)){
          # toxicity high
          if(prob_t_sucess[[i]][xt]<1-int_toxicity_prob){
            if(int_toxicity_prob!=1){
            decision[[i]][xr,xt]=1
            }
            # Futility high
          } else if(prob_r_sucess[[i]][xr]<1-int_futility_prob){
            decision[[i]][xr,xt]=1
            # combined high
          } else if(prob_r_sucess[[i]][xr]*prob_t_sucess[[i]][xt]<1-int_combined_prob){
            decision[[i]][xr,xt]=1
            #efficacy and toxicity good
          } else if(prob_r_sucess[[i]][xr]*prob_t_sucess[[i]][xt]>int_efficacy_prob){
            decision[[i]][xr,xt]=3
          } else {
            decision[[i]][xr,xt]=2
          }
        }
      }
    }
  }

  # final analysis
  decision[[length(reviews)]]=matrix(1,max(reviews)+1,max(reviews)+1)
  for(xr in 1:(max(reviews)+1)){
    for(xt in 1:(max(reviews)+1)){
      if((xr-1)>=R_cut & (xt-1)<=T_cut){
      decision[[length(reviews)]][xr,xt]=3
      }
    }
  }


  # Cutpoints update from decision matrices
  for(i in 1:length(reviews)){
    ans=which(apply(decision[[i]]==3,2,sum)>0)
    cutpoints[i,2]=ifelse(length(ans)>0,max(ans)-1,NA)
    cutpoints[i,2]=ifelse(cutpoints[i,2]<0,NA,cutpoints[i,2])
    ans=which(apply(decision[[i]]==2,2,sum)>0)
    cutpoints[i,3]=ifelse(length(ans)>0,max(ans),NA)
    ans=which(apply(decision[[i]]==2,1,sum)>0)
    cutpoints[i,4]=ifelse(length(ans),min(ans)-2,NA)
    cutpoints[i,4]=ifelse(cutpoints[i,4]<0,NA,cutpoints[i,4])
    ans=which(apply(decision[[i]]==3,1,sum)>0)
    cutpoints[i,5]=ifelse(length(ans)>0,min(ans)-1,NA)
  }

  # stop rules for final analysis
  cutpoints[i,3]=cutpoints[i,2]+1
  cutpoints[i,4]=cutpoints[i,5]-1

  # Frequentist analysis of decision
  # Prepare data for return and print
  # Bayesian analysis of decision
  return(.properties_binom_two(t,r,reviews,pra,prb,pta,ptb,futility_critical_value,efficacy_critical_value,toxicity_critical_value,no_toxicity_critical_value,cutpoints,decision))

}

# all done
