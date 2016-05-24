# LINES trial exact "simulation" of expected outcome under specific data conditions.
# Since the number of patients is small a probabilistic approach is more efficient than a large sampling approach.
# Loss function based methodology
# This is based on the paper by Yiyi Chen 2009 Adaptive group sequential design for phase II clinical trials;
# This version does not trade-off efficacy and toxicity and considers a square section for superiority.

# if stage_after_trial is 0 or missing you force the trial to decide between sucess or failure based on the final number of patients.



bayes_binom_two_loss=function(t,r,reviews,pra=1,prb=1,pta=1,ptb=1,l_alpha_beta,l_alpha_c,stage_after_trial=0,fun.integrate,efficacy_critical_value=0.2,toxicity_critical_value=0.1,futility_critical_value=0.35,no_toxicity_critical_value=0.3,decision=NULL,W=NULL,fun.graph=NULL,...){

  if(length(l_alpha_c)!=max(reviews,stage_after_trial)){
    print("The cost function is constant for all patients")
    l_alpha_c=rep(l_alpha_c[1],max(reviews,stage_after_trial))
  }
  l_beta=1
  l_alpha=l_alpha_beta*l_beta
  c=l_alpha/l_alpha_c


  # variable to store the critical values for each endpoint
  cutpoints=matrix(NA,length(reviews),5)
  cutpoints[,1]=reviews

  if(max(stage_after_trial)>max(reviews)){
    stages=c(reviews,stage_after_trial)
  } else {
    stages=reviews
  }


  ##################################################################################
  # Compute all decision rules and probabilities
  ##################################################################################
  # given prior,xr,xt and n pi can be calculated

  # if decision rule and W already provided skip this section
  if(is.null(decision) | is.null(decision)){
    W=list()

    for(i in 1:length(stages)){
      W[[i]]=matrix(0,stages[i]+1,stages[i]+1)
      for(xr in 0:stages[i]){
        for(xt in 0:stages[i]){
          # posterior distribution
          ar=pra+xr
          br=prb+stages[i]-xr
          at=pra+xt
          bt=prb+stages[i]-xt
          W[[i]][xr+1,xt+1]=fun.integrate(ar,br,at,bt,...)
        }
      }
    }

    # calculate the decision rule for each xr,xt and review
    exp_loss_d3_k1=list()
    decision=list()
    for(i in length(stages):1){

      exp_loss_d3_k1[[i]]=matrix(0,stages[i]+1,stages[i]+1)
      decision[[i]]=matrix(0,stages[i]+1,stages[i]+1)
      for(xr in 0:stages[i]){
        for(xt in 0:stages[i]){


          # expected futility loss function for stage K
          exp_loss_d1_k=sum(c[1:stages[i]])+W[[i]][xr+1,xt+1]*l_beta
          # expected efficacy loss function for stage K
          exp_loss_d2_k=sum(c[1:stages[i]])+(1-W[[i]][xr+1,xt+1])*l_alpha

          # future outcomes
          # compute the probability of each pair of xr_k+1,xt_k+1 occuring using posterior predictive distribution
          # response part of posterior predictive
          if((i+1)>length(stages)){
            # make a definative decision between d1 and d2 not comparing to a future analysis

            exp_loss_d3_k1[[i]][xr+1,xt+1]=min(exp_loss_d1_k,exp_loss_d2_k)
            decision[[i]][xr+1,xt+1]=3-2*(exp_loss_d1_k<exp_loss_d2_k)
          } else if((i+1)==length(stages)){
            # final stage but allowed to recommend further analysis to a future analysis
            a3=rep(0,stages[i+1]+1)
            a3[(xr+1):(xr+1+stages[i+1]-stages[i])]=dbetabinom.ab(0:(stages[i+1]-stages[i]),stages[i+1]-stages[i],pra+xr,prb+stages[i]-xr)
            # toxicity part of posterior predictive
            a4=rep(0,stages[i+1]+1)
            a4[(xt+1):(xt+1+stages[i+1]-stages[i])]=dbetabinom.ab(0:(stages[i+1]-stages[i]),stages[i+1]-stages[i],pta+xt,ptb+stages[i]-xt)
            post_pred=a3%*%t(a4)

            exp_loss_d1_k1=sum(c[1:stages[i+1]])+W[[i+1]]*l_beta
            exp_loss_d2_k1=sum(c[1:stages[i+1]])+(1-W[[i+1]])*l_alpha
            exp_loss_d3_k=sum(post_pred*pmin(exp_loss_d1_k1,exp_loss_d2_k1))




            exp_loss_d3_k1[[i]][xr+1,xt+1]=min(exp_loss_d1_k,exp_loss_d2_k,exp_loss_d3_k)
            decision[[i]][xr+1,xt+1]=3-1*(exp_loss_d3_k<min(exp_loss_d1_k,exp_loss_d2_k))-2*(exp_loss_d3_k>min(exp_loss_d1_k,exp_loss_d2_k))*(exp_loss_d1_k<exp_loss_d2_k)
          } else {
            # interim stage combine future analysis to this analysis
            a3=rep(0,stages[i+1]+1)
            a3[(xr+1):(xr+1+stages[i+1]-stages[i])]=dbetabinom.ab(0:(stages[i+1]-stages[i]),stages[i+1]-stages[i],pra+xr,prb+stages[i]-xr)
            # toxicity part of posterior predictive
            a4=rep(0,stages[i+1]+1)
            a4[(xt+1):(xt+1+stages[i+1]-stages[i])]=dbetabinom.ab(0:(stages[i+1]-stages[i]),stages[i+1]-stages[i],pta+xt,ptb+stages[i]-xt)
            post_pred=a3%*%t(a4)
            exp_loss_d1_k1=sum(c[1:stages[i+1]])+W[[i+1]]*l_beta
            exp_loss_d2_k1=sum(c[1:stages[i+1]])+(1-W[[i+1]])*l_alpha
            #exp_loss_d3_k1[[i]] from previous loop at k+1 level
            exp_loss_d3_k=sum(post_pred*pmin(exp_loss_d1_k1,exp_loss_d2_k1,exp_loss_d3_k1[[i+1]]))

            exp_loss_d3_k1[[i]][xr+1,xt+1]=min(exp_loss_d1_k,exp_loss_d2_k,exp_loss_d3_k)
            decision[[i]][xr+1,xt+1]=3-1*(exp_loss_d3_k<min(exp_loss_d1_k,exp_loss_d2_k))-2*(exp_loss_d3_k>min(exp_loss_d1_k,exp_loss_d2_k))*(exp_loss_d1_k<exp_loss_d2_k)
          }
        }
      }
    }
  }
  ##################################################################################


  # Cutpoints update from decision matrices
  for(i in 1:length(reviews)){
    ans=which(apply(decision[[i]]==3,2,sum)>0)
    cutpoints[i,2]=ifelse(length(ans)>0,max(ans)-1,NA)
    cutpoints[i,2]=ifelse(cutpoints[i,2]<0,NA,cutpoints[i,2])
    ans=which(apply(decision[[i]] == 2 | decision[[i]] == 3,2,sum)>0)
    cutpoints[i,3]=ifelse(length(ans)>0,max(ans),NA)
    ans=which(apply(decision[[i]] == 2 | decision[[i]] == 3,1,sum)>0)
    cutpoints[i,4]=ifelse(length(ans),min(ans)-2,NA)
    cutpoints[i,4]=ifelse(cutpoints[i,4]<0,NA,cutpoints[i,4])
    ans=which(apply(decision[[i]]==3,1,sum)>0)
    cutpoints[i,5]=ifelse(length(ans)>0,min(ans)-1,NA)
  }


  # Frequentist analysis of decision
  # Prepare data for return and print
  # Bayesian analysis of decision
  to.return=.properties_binom_two(t,r,reviews,pra,prb,pta,ptb,futility_critical_value,efficacy_critical_value,toxicity_critical_value,no_toxicity_critical_value,cutpoints,decision)

  to.return@graph=list(fun.graph=fun.graph,max.patients=max(stage_after_trial,reviews),...)

  return(to.return)

}





#####################################################################################
#####################################################################################
## Integrate over the region
#####################################################################################
#####################################################################################


################################################################################
# intercept of ratio line and square region
tradeoff_ratio_intercepts=function(R_min,R_max,T_min,T_max,theta=0){
  qt=(T_min+T_max)/2
  qr=(R_min+R_max)/2
  pt=T_min
  pr=R_min
  if(theta==0){theta=(qr*(1-qt))/((1-qr)*qt)}

  prmin=max(pt*theta/(1-pt+pt*theta),R_min)
  ptmin=max(pr/(theta-pr*theta+pr),T_min)
  pt=T_max
  pr=R_max
  prmax=min(pt*theta/(1-pt+pt*theta),R_max)
  ptmax=min(pr/(theta-pr*theta+pr),T_max)
  return(list(prmin=prmin,prmax=prmax,ptmin=ptmin,ptmax=ptmax))
}

# integrate over acceptable region
# Odds ratio region
tradeoff_ratio_integrate=function(ar,br,at,bt,efficacy_region_min,toxicity_region_max,efficacy_region_max,toxicity_region_min,theta,intercepts){
  # Function to integrate over tradeoff region

  vec.integrate.over.trade=function(p,ar,br,at,bt,theta){
    res=rep(0,length(p))
    for(i in 1:length(p)){
      res[i]=(1-pbeta(p[i]*theta/(1-theta*p[i]+p[i]),ar,br))*dbeta(p[i],at,bt)
      #
    }
    return(res)
  }

  return(((1-pbeta(intercepts$prmin,ar,br))*pbeta(intercepts$ptmin,at,bt))+(integrate(vec.integrate.over.trade,intercepts$ptmin,intercepts$ptmax,ar,br,at,bt,theta)$value)+((1-pbeta(intercepts$prmax,ar,br))*(pbeta(toxicity_region_max,at,bt)-pbeta(intercepts$ptmax,at,bt))))

}



################################################################################
# integrate over acceptable region
# square region
tradeoff_square_integrate=function(ar,br,at,bt,efficacy_region_min,toxicity_region_max){
  return((1-pbeta(efficacy_region_min,ar,br))*pbeta(toxicity_region_max,at,bt))
}


################################################################################
# integrate over acceptable region
# ellipse region
tradeoff_ellipse_integrate=function(ar,br,at,bt,efficacy_region_min,toxicity_region_max,efficacy_region_max,toxicity_region_min){
  # Function to integrate over tradeoff region
  vec.integrate.over.trade=function(p,ar,br,at,bt){
    res=rep(0,length(p))
    for(i in 1:length(p)){
      res[i]=(1-pbeta(-sqrt(1-((p[i]-toxicity_region_min)/(toxicity_region_max-toxicity_region_min))^2)*(efficacy_region_max-efficacy_region_min)+efficacy_region_max,ar,br))*dbeta(p[i],at,bt)
      # pr=-sqrt(1-((pt-toxicity_region_min)/(toxicity_region_max-toxicity_region_min))^2)*(efficacy_region_max-efficacy_region_min)+efficacy_region_max
    }
    return(res)
  }

  return(((1-pbeta(efficacy_region_min,ar,br))*pbeta(toxicity_region_min,at,bt))+integrate(vec.integrate.over.trade,toxicity_region_min,toxicity_region_max,ar,br,at,bt)$value)

}



################################################################################
# integrate over acceptable region
# Linear tradeoff
tradeoff_linear_integrate=function(ar,br,at,bt,efficacy_region_min,toxicity_region_max,efficacy_region_max,toxicity_region_min){
  vec.integrate.over.trade=function(p,ar,br,at,bt){
    res=rep(0,length(p))
    for(i in 1:length(p)){
      res[i]=(1-pbeta((efficacy_region_max-efficacy_region_min)/(toxicity_region_max-toxicity_region_min)*p[i]+efficacy_region_min-((efficacy_region_max-efficacy_region_min)/(toxicity_region_max-toxicity_region_min))*toxicity_region_min,ar,br))*dbeta(p[i],at,bt)
    }
    return(res)
  }

  return(((1-pbeta(efficacy_region_min,ar,br))*pbeta(toxicity_region_min,at,bt))+integrate(vec.integrate.over.trade,toxicity_region_min,toxicity_region_max,ar,br,at,bt)$value)

}


#####################################################################################
#####################################################################################
## Plot graph of region
#####################################################################################
#####################################################################################

tradeoff_ratio_graph=function(input){

  ##################################################################################
  # Graph of acceptable region
  pt=0:1000/1000
  pr=pt*input$theta / (1-pt+input$theta*pt )
  plot(0,0,xlim=c(0,1),ylim=c(0,1),col="white",xlab="response",ylab="toxicity", yaxs="i",xaxs="i",main=paste("Maximum for backward induction =" ,input$max.patients))

  polygon(c(input$efficacy_region_min,1,1,input$efficacy_region_min),c(input$toxicity_region_min,input$toxicity_region_min,0,0),col="gray")
  polygon(c(input$efficacy_region_max,1,1,input$efficacy_region_max),c(input$toxicity_region_min,input$toxicity_region_min,input$toxicity_region_max ,input$toxicity_region_max ),col="gray")
  polygon(c(input$efficacy_region_min,pr[which(pt>=input$toxicity_region_min & pt<=input$toxicity_region_max  & pr>=input$efficacy_region_min & pr<=input$efficacy_region_max)],input$efficacy_region_max,input$efficacy_region_max),c(input$toxicity_region_min,pt[which(pt>=input$toxicity_region_min & pt<=input$toxicity_region_max  & pr>=input$efficacy_region_min & pr<=input$efficacy_region_max)],input$toxicity_region_max ,input$toxicity_region_min),col="gray")

  abline(h=c(input$toxicity_region_min,input$toxicity_region_max))
  abline(v=c(input$efficacy_region_min,input$efficacy_region_max))
  lines(pr,pt)

}


tradeoff_square_graph=function(input){
  ##################################################################################
  # Graph of acceptable region for square region
  plot(0,0,xlim=c(0,1),ylim=c(0,1),col="white",xlab="response",ylab="toxicity", yaxs="i",xaxs="i",main=paste("Maximum for backward induction =" ,input$max.patients))
  polygon(c(input$efficacy_region_min,1,1,input$efficacy_region_min),c(input$toxicity_region_max,input$toxicity_region_max,0,0),col="gray")
  abline(h=input$toxicity_region_max)
  abline(v=input$efficacy_region_min)
}


tradeoff_ellipse_graph=function(input){
  ##################################################################################
  # Graph of acceptable region
  pr=((input$efficacy_region_min*1000):(input$efficacy_region_max*1000))/1000
  pt=sqrt(1-((pr-input$efficacy_region_max)/(input$efficacy_region_max-input$efficacy_region_min))^2)*(input$toxicity_region_max-input$toxicity_region_min)+input$toxicity_region_min
  plot(0,0,xlim=c(0,1),ylim=c(0,1),col="white",xlab="Response",ylab="Toxicity", yaxs="i",xaxs="i",main=paste("Maximum for backward induction =" ,input$max.patients))
  polygon(c(pr,1,1,input$efficacy_region_min),c(pt,input$toxicity_region_max,0,0),col="gray")
  abline(h=input$toxicity_region_min)
  abline(h=input$toxicity_region_max)
  abline(v=input$efficacy_region_min)
  abline(v=input$efficacy_region_max)
  lines(pr,pt)

}



tradeoff_linear_graph=function(input){
  ##################################################################################
  # Graph of acceptable region for square region
  plot(0,0,xlim=c(0,1),ylim=c(0,1),col="white",xlab="response",ylab="toxicity", yaxs="i",xaxs="i",main=paste("Maximum for backward induction =" ,input$max.patients))
  polygon(c(input$efficacy_region_min,input$efficacy_region_max,1,1,input$efficacy_region_min),c(input$toxicity_region_min,input$toxicity_region_max,input$toxicity_region_max,0,0),col="gray")
  abline(h=input$toxicity_region_min)
  abline(h=input$toxicity_region_max)
  abline(v=input$efficacy_region_min)
  abline(v=input$efficacy_region_max)
}
