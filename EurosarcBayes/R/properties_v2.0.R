##################################################################


properties=function(x,...){

  cat("\nThere is no default function for this function type.")
  cat("\nThe object must have a class attached to reach the required function")
}


setGeneric("properties")

##################################################################
##################################################################
##################################################################


# Generate trial properties for a binomial single arm design
properties_binom_one=function(failure=NULL,success=NULL,reviews=NULL,p0,p1,prior.a=0,prior.b=0,round=TRUE,cutpoints=NULL){

  if(!is.null(cutpoints)){
    reviews=cutpoints$n
    success=cutpoints$efficacy
    failure=cutpoints$futility
  }

  alpha=0
  beta=0
  ln=length(reviews)
  maxn=max(reviews)

  # if missing generate
  if(is.null(failure)){
    failure=rep(NA,length(success))
  }
  if(is.null(success)){
    success=rep(NA,length(failure))
  }

  # if last stop for futility not given assign it for easy calculation
  if(is.na(failure[length(failure)])){
    failure[length(failure)]=success[length(success)]-1
  }
  # if last stop for efficacy not given assign it for easy calculation
  if(is.na(success[length(success)])){
    failure[length(success)]=failure[length(success)]+1
  }

  success=as.numeric(success)
  failure=as.numeric(failure)

  ###################################################################
  # frequentist properties
  alpha=0
  beta=0
  exp.p0=0
  exp.p1=0
  old0=c(1,rep(0,max(reviews)))
  old1=c(1,rep(0,max(reviews)))
  for(i in 1:ln){
    if(i>1){
      patients=reviews[i-1]
    } else {
      patients=0
    }

    new0=rep(0,max(reviews)+1)
    new1=rep(0,max(reviews)+1)

    # draw patients from distribution to get to stage i.
    probp0=dbinom(0:(reviews[i]-patients),reviews[i]-patients,p0)
    probp1=dbinom(0:(reviews[i]-patients),reviews[i]-patients,p1)
    # update current trial distribution
    for(j in 1:(patients+1)){
      new0[j:(j+reviews[i]-patients)]=new0[j:(j+reviews[i]-patients)] + old0[j]*probp0
      new1[j:(j+reviews[i]-patients)]=new1[j:(j+reviews[i]-patients)] + old1[j]*probp1
    }


    stop.0=0
    stop.1=0
    if(!is.null(failure)){
      if(!is.na(failure[i])){
        stop.0=stop.0+sum(new0[1:(failure[i]+1)])
        stop.1=stop.1+sum(new1[1:(failure[i]+1)])
        beta=beta+sum(new1[1:(failure[i]+1)])
        new0[1:(failure[i]+1)]=0
        new1[1:(failure[i]+1)]=0
      }}

    if(!is.na(success[i])){
      stop.0=stop.0+sum(new0[(success[i]+1):(reviews[i]+1)])
      stop.1=stop.1+sum(new1[(success[i]+1):(reviews[i]+1)])
      alpha=alpha+sum(new0[(success[i]+1):(reviews[i]+1)])
      new0[(success[i]+1):(reviews[i]+1)]=0
      new1[(success[i]+1):(reviews[i]+1)]=0
    }

    exp.p0=exp.p0+reviews[i]*stop.0
    exp.p1=exp.p1+reviews[i]*stop.1
    old0=new0
    old1=new1
  }

  ###################################################################
  # Bayesian properties
  eta=pbeta(p1,prior.a+failure,prior.b+reviews-failure)
  zeta=1-pbeta(p0,prior.a+success,prior.b+reviews-success)


  if(round==TRUE){
    alpha=round(alpha,3)
    beta=round(beta,3)
    exp.p0=round(exp.p0,2)
    exp.p1=round(exp.p1,2)
    eta=round(eta,3)
    zeta=round(zeta,3)
  }

  return(trialDesign_binom_one(p0=p0,p1=p1,reviews=reviews,failure=failure,success=success,alpha=alpha,power=1-beta,exp.p0=exp.p0,exp.p1=exp.p1,eta=eta,zeta=zeta))

}

##################################################################
##################################################################
##################################################################
#output$design
setMethod(f="properties",signature="binom_two_bryantday",definition=function(x,t,r,pra,prb,pta,ptb,futility_critical_value=0.2,efficacy_critical_value=0.35,toxicity_critical_value=0.3,no_toxicity_critical_value=0.1,target="optimal"){

  if(target=="optimal"){
    output=x@optimal
  } else if(target=="minmax"){
    output=x@minmax
  } else {
    output=x@all.fit[,target]
  }

  reviews=c(output$Interim,output$final)

  # create cutpoint matrix
  cutpoints=matrix(NA,2,5)
  cutpoints[,1]=c(output$Interim,output$final)
  cutpoints[2,2]=cutpoints[2,1]-output$toxicity
  cutpoints[1,3]=cutpoints[1,1]-output$Itoxicity+1
  cutpoints[2,3]=cutpoints[2,1]-output$toxicity+1
  cutpoints[1,4]=output$Iresponse-1
  cutpoints[2,4]=output$response-1
  cutpoints[2,5]=output$response


  # Frequentist analysis of decision
  # Prepare data for return and print
  # Bayesian analysis of decision
  return(.properties_binom_two(t,r,reviews,pra,prb,pta,ptb,futility_critical_value,efficacy_critical_value,toxicity_critical_value,no_toxicity_critical_value,cutpoints,decision=NULL))

}
)

##################################################################
##################################################################
##################################################################
setMethod(f="properties",signature="binom_two_singlestage",definition=function(x,t,r,pra,prb,pta,ptb,futility_critical_value=0.2,efficacy_critical_value=0.35,toxicity_critical_value=0.3,no_toxicity_critical_value=0.1){


  output=x@optimal
  reviews=output$n
  cutpoints=matrix(NA,1,5)
  cutpoints[1,1]=output$n
  cutpoints[1,2]=output$toxicity
  cutpoints[1,3]=output$toxicity+1
  cutpoints[1,4]=output$response-1
  cutpoints[1,5]=output$response

  # Frequentist analysis of decision
  # Prepare data for return and print
  # Bayesian analysis of decision

  return(.properties_binom_two(t,r,reviews,pra,prb,pta,ptb,futility_critical_value,efficacy_critical_value,toxicity_critical_value,no_toxicity_critical_value,cutpoints,decision=NULL))

}
)

# ended
##################################################################
##################################################################
##################################################################


.properties_binom_two=function(t,r,reviews,pra,prb,pta,ptb,futility_critical_value,efficacy_critical_value,toxicity_critical_value,no_toxicity_critical_value,cutpoints,decision=NULL){

  if(is.null(decision)){
    le=function(x,y){
      if(is.na(y)){
        return(FALSE)
      } else {
        return(x<=y)
      }
    }

    ge=function(x,y){
      if(is.na(y)){
        return(FALSE)
      } else {
        return(x>=y)
      }
    }

    ########################################################################
    # create decision matrices for comparison to other methods
    decision=list()
    for(i in 1:length(reviews)){
      decision[[i]]=matrix(2,reviews[i]+1,reviews[i]+1)
      for(xr in 0:(reviews[i])){
        for(xt in 0:(reviews[i])){
          if(ge(xr,cutpoints[i,5]) & le(xt,cutpoints[i,2])){
            decision[[i]][xr+1,xt+1]=3
          } else if(ge(xt,cutpoints[i,3]) | le(xr,cutpoints[i,4])){
            decision[[i]][xr+1,xt+1]=1
          }
        }
      }
    }
  }



  no_reviews=length(reviews)
  # Lines Bayesian methods.
  # Run frequentist and bayesian analysis of design.
  # prepare for return
  ########################################

  # prepare variable to confirm probability is preserved
  precision=rep(0,length(r))

  table=matrix(,6,length(r)+1)
  table[1,1]="Stop early - Futility/Toxicity"
  table[2,1]="Stop early - Efficacy"
  table[3,1]="Continue to final analysis - Efficacy"
  table[4,1]="Continue to final analysis - Futility/Toxicity"
  table[5,1]="Continue to final analysis, recommend a further stage (inconclusive)"
  table[6,1]="Expected number of patients recruited"


  results=list()
  # Run for each different scenario
  for(m in 1:length(r)){

    # set up outcome matrix with initial state 0 patients.
    outcomer=rep(0,max(reviews)+1)
    outcomet=rep(0,max(reviews)+1)
    outcomer[1]=1
    outcomet[1]=1

    # prepare tempoury variables for outcome
    noutcomer=rep(0,max(reviews)+1)
    noutcomet=rep(0,max(reviews)+1)

    # variable to save the results for early exit
    results[[m]]=matrix(0,0,3)

    patients=0
    for(i in 1:no_reviews){
      #print(c(m,reviews[i]))

      # simulate from operating characteristics
      for(j in 1:(patients+1)){
        noutcomer[j:(j+reviews[i]-patients)]=noutcomer[j:(j+reviews[i]-patients)]+outcomer[j]*dbinom(0:(reviews[i]-patients),reviews[i]-patients,r[m])
        noutcomet[j:(j+reviews[i]-patients)]=noutcomet[j:(j+reviews[i]-patients)]+outcomet[j]*dbinom(0:(reviews[i]-patients),reviews[i]-patients,t[m])
      }


      total=sum(noutcomer)
      outcomer=rep(0,max(reviews)+1)
      outcomet=rep(0,max(reviews)+1)

      # record the results and stop those trials which stop early.
      #two lines for early stopping
      results[[m]]=rbind(results[[m]],0)
      results[[m]]=rbind(results[[m]],0)
      if(reviews[i]==max(reviews)){
        results[[m]]=rbind(results[[m]],0)
      }
      for(xr in 0:reviews[i]){
        # skip if probability is zero
        if(noutcomer[xr+1]!=0){
          for(xt in 0:reviews[i]){
            # skip if probability is zero
            if(noutcomet[xt+1]!=0){

              # If not the final review
              if(reviews[i]!=max(reviews)){
                if(decision[[i]][xr+1,xt+1]==1){
                  results[[m]][dim(results[[m]])[1]-1,]=c(reviews[i],1,results[[m]][dim(results[[m]])[1]-1,3]+noutcomer[xr+1]*noutcomet[xt+1]/total)
                  # stop for futility
                } else if(decision[[i]][xr+1,xt+1]==3){
                  results[[m]][dim(results[[m]])[1],]=c(reviews[i],2,results[[m]][dim(results[[m]])[1],3]+noutcomer[xr+1]*noutcomet[xt+1]/total)
                  # stop for efficacy
                } else if(decision[[i]][xr+1,xt+1]==2){
                  outcomer[xr+1]=outcomer[xr+1]+noutcomer[xr+1]*noutcomet[xt+1]/total
                  outcomet[xt+1]=outcomet[xt+1]+noutcomer[xr+1]*noutcomet[xt+1]/total
                  # continue to next stage
                }

                # If the final review
              } else {
                if(decision[[i]][xr+1,xt+1]==1){
                  results[[m]][dim(results[[m]])[1]-2,]=c(reviews[i],4,results[[m]][dim(results[[m]])[1]-2,3]+noutcomer[xr+1]*noutcomet[xt+1]/total)
                  # stop for futility
                } else if(decision[[i]][xr+1,xt+1]==3){
                  results[[m]][dim(results[[m]])[1]-1,]=c(reviews[i],3,results[[m]][dim(results[[m]])[1]-1,3]+noutcomer[xr+1]*noutcomet[xt+1]/total)
                  # stop for efficacy
                } else if(decision[[i]][xr+1,xt+1]==2){
                  results[[m]][dim(results[[m]])[1],]=c(reviews[i],5,results[[m]][dim(results[[m]])[1],3]+noutcomer[xr+1]*noutcomet[xt+1]/total)
                  # continue to next stage
                  # But this is the final allowed stage
                }
              }
              #print(c(stages[i],decision,exp_loss_d3_k,exp_loss_d4_k))
              #############################################################

            }
          }
        }
      }

      noutcomer=rep(0,max(reviews)+1)
      noutcomet=rep(0,max(reviews)+1)
      patients=reviews[i]
    }

    # store results to the table
    table[1,m+1]=sum(subset(results[[m]][,3],results[[m]][,2]==1))*100
    table[2,m+1]=sum(subset(results[[m]][,3],results[[m]][,2]==2))*100
    table[3,m+1]=sum(subset(results[[m]][,3],results[[m]][,2]==3))*100
    table[4,m+1]=sum(subset(results[[m]][,3],results[[m]][,2]==4))*100
    table[5,m+1]=sum(subset(results[[m]][,3],results[[m]][,2]==5))*100
    table[6,m+1]=sum(results[[m]][,1]*results[[m]][,3])

    # total sum (any difference is an R precision problem (provided coding is correct)
    precision[m]=sum(results[[m]][,3])

  }

  ################################################################################
  ################################################################################
  ################################################################################


  # Lines Bayesian methods.
  # Preparation for return.
  ################################################################################
  # Prepare tables for output
  data=data.frame("x"=table[,1])
  st=c("Stopping rules")
  for(i in 2:dim(table)[2]){
    data=data.frame(data,x=round(as.numeric(table[,i]),2))
    st=c(st,paste0("T=",t[i-1],", R=",r[i-1]))
  }

  data=setnames(data,st)

  cutpoints=data.frame(cutpoints)
  cutpoints=rename(cutpoints,c("X1"="patient review", "X2"="low toxicity", "X3"="high toxicity", "X4"="poor outcome", "X5"="good outcome"))

  cat("cut-points at each analysis\n")
  print(cutpoints)
  cat("\nFrequentist properties of design\n")
  data=data[c(1,4,2,3,6),]
  print(data)


  if(is.na(cutpoints[no_reviews,3])){cutpoints[no_reviews,3]=cutpoints[no_reviews,2]+1}
  if(is.na(cutpoints[no_reviews,4])){cutpoints[no_reviews,4]=cutpoints[no_reviews,5]-1}

  bayes.properties=matrix(0,nrow=no_reviews,ncol=9)

  for(i in 1:no_reviews){
    bayes.properties[i,]=c(cutpoints[i,1],1-pbeta(no_toxicity_critical_value,pta+cutpoints[i,2],ptb+cutpoints[i,1]-cutpoints[i,2]),1-pbeta(toxicity_critical_value,pta+cutpoints[i,2],ptb+cutpoints[i,1]-cutpoints[i,2]),
                           1-pbeta(no_toxicity_critical_value,pta+cutpoints[i,3],ptb+cutpoints[i,1]-cutpoints[i,3]),1-pbeta(toxicity_critical_value,pta+cutpoints[i,3],ptb+cutpoints[i,1]-cutpoints[i,3]),
                           1-pbeta(efficacy_critical_value,pra+cutpoints[i,4],prb+cutpoints[i,1]-cutpoints[i,4]),1-pbeta(futility_critical_value,pra+cutpoints[i,4],prb+cutpoints[i,1]-cutpoints[i,4]),
                           1-pbeta(efficacy_critical_value,pra+cutpoints[i,5],prb+cutpoints[i,1]-cutpoints[i,5]),1-pbeta(futility_critical_value,pra+cutpoints[i,5],prb+cutpoints[i,1]-cutpoints[i,5]))
  }



  bayes.properties.dta=data.frame(matrix(0,length(reviews),9))
  bayes.properties.dta[,1]=bayes.properties[,1]
  for(i in 2:9){
    bayes.properties.dta[,i]=round(bayes.properties[,i],3)
  }


  bayes.properties.dta=setnames(bayes.properties.dta,c("n",paste0("T>",no_toxicity_critical_value),paste0("T>",toxicity_critical_value),paste0("T>",no_toxicity_critical_value),paste0("T>",toxicity_critical_value),paste0("R>",efficacy_critical_value),paste0("R>",futility_critical_value),paste0("R>",efficacy_critical_value),paste0("R>",futility_critical_value)))

  cat("\nBayesian properties of trial design\n")
  print(bayes.properties.dta,row.names = FALSE)

  futility=1-max(bayes.properties.dta[,7],na.rm = T)
  efficacy=min(bayes.properties.dta[,8],na.rm = T)
  no_toxicity=ifelse(sum(is.na(bayes.properties.dta[,2]))==length(bayes.properties.dta[,2]),NA,1-max(bayes.properties.dta[,2],na.rm = T))
  toxicity=ifelse(sum(is.na(bayes.properties.dta[,5]))==length(bayes.properties.dta[,5]),NA,min(bayes.properties.dta[,5],na.rm = T))
  cat(paste0("\nFutility     P(R<",futility_critical_value,")=",futility))
  cat(paste0("\nEfficacy     P(R>",efficacy_critical_value,")=",efficacy,"\n"))
  cat(paste0("\nToxicity ok  P(T<",no_toxicity_critical_value,")=",no_toxicity))
  cat(paste0("\nToxicity     P(T>",toxicity_critical_value,")=",toxicity))

  return(trialDesign_binom_two(reviews=reviews,data=data,cutpoints=cutpoints,precision=precision,decision=decision,post.futility=futility, post.efficacy=efficacy, post.toxicity=toxicity,post.no_toxicity=no_toxicity,graph=list()))

}
