# LINES trial exact "simulation" of expected outcome under specific data conditions.
# Since the number of patients is small a probabilistic approach is more efficient than a large sampling approach.
# play with stopping rules in the way of the first trial.





bayes_binom_two_postprob=function(t,r,reviews,pra=1,prb=1,pta=1,ptb=1,futility_critical_value,futility_prob_stop,efficacy_critical_value,efficacy_prob_stop,toxicity_critical_value,toxicity_prob_stop,no_toxicity_critical_value,no_toxicity_prob_stop=NULL){

  no_reviews=length(reviews)
  if(length(futility_prob_stop)==1){
    futility_prob_stop=rep(futility_prob_stop,no_reviews)
  } else if(length(futility_prob_stop)!=no_reviews){
    stop("The length of futility_prob_stop does not match number of reviews")
  }
  if(length(efficacy_prob_stop)==1){
    efficacy_prob_stop=rep(efficacy_prob_stop,no_reviews)
  } else if(length(efficacy_prob_stop)!=no_reviews){
    stop("The length of efficacy_prob_stop does not match number of reviews")
  }
  if(length(toxicity_prob_stop)==1){
    toxicity_prob_stop=rep(toxicity_prob_stop,no_reviews)
  } else if(length(toxicity_prob_stop)!=no_reviews){
    stop("The length of toxicity_prob_stop does not match number of reviews")
  }
  if(length(no_toxicity_prob_stop)==1){
    no_toxicity_prob_stop=rep(no_toxicity_prob_stop,no_reviews)
  } else if(length(no_toxicity_prob_stop)!=no_reviews){
    stop("The length of no_toxicity_prob_stop does not match number of reviews")
  }
  if(length(pra)!=1 | length(prb)!=1 | length(pta)!=1 | length(ptb)!=1){
    stop("The length of the prior parameters is not one each")
  }

  # variable to store the critical values for each endpoint
  cutpoints=matrix(NA,no_reviews,5)

  ########################################################################
  # compute cutpoints. THis does all the Bayesian work to design this method!
  for(i in 1:no_reviews){

    vals=which(pbeta(toxicity_critical_value,pta+0:reviews[i],ptb+reviews[i]-0:reviews[i])<1-toxicity_prob_stop[i])
    hightox=ifelse(length(vals)>0,min(vals)-1,NA)
    vals=which(pbeta(futility_critical_value,pra+0:reviews[i],prb+reviews[i]-0:reviews[i])>futility_prob_stop[i])
    pooroutcome=ifelse(length(vals)>0,max(vals)-1,NA)
    vals=which(pbeta(efficacy_critical_value,pra+0:reviews[i],prb+reviews[i]-0:reviews[i])<1-efficacy_prob_stop[i])
    goodoutcome=ifelse(length(vals)>0,min(vals)-1,NA)
    vals=which(pbeta(no_toxicity_critical_value,pta+0:reviews[i],ptb+reviews[i]-0:reviews[i])>no_toxicity_prob_stop[i])
    lowtox=ifelse(length(vals)>0,max(vals)-1,NA)
   if(!is.na(lowtox) & !is.na(hightox)){
     if(lowtox>=hightox){
       hightox=lowtox+1
     }
    }
    if(!is.na(pooroutcome) & !is.na(goodoutcome)){
      if(goodoutcome<=pooroutcome){
        pooroutcome=goodoutcome-1
      }
    }
   if(is.na(lowtox) & !is.na(goodoutcome)){
     lowtox=hightox-1
   }
    cutpoints[i,]=c(reviews[i],lowtox,hightox,pooroutcome,goodoutcome)
  }
  if(is.na(cutpoints[no_reviews,3])){
    cutpoints[no_reviews,3]=cutpoints[no_reviews,2]+1
  }
  if(is.na(cutpoints[no_reviews,4])){
    cutpoints[no_reviews,4]=cutpoints[no_reviews,5]-1
  }



  # Frequentist analysis of decision
  # Prepare data for return and print
  # Bayesian analysis of decision
  return(.properties_binom_two(t,r,reviews,pra,prb,pta,ptb,futility_critical_value,efficacy_critical_value,toxicity_critical_value,no_toxicity_critical_value,cutpoints,decision=NULL))

}

# all done


