
freq_binom_one_landemets=function(reviews,p0,p1,r=c(p0,p1),alpha=0.1,beta=0.1,prior.a=0,prior.b=0){

  t=reviews/reviews[length(reviews)]

  # Obrien flemming bounds
  probcuta=2*(1-pnorm((qnorm(1-alpha/2))/sqrt(t)))
  probcutb=2*(1-pnorm((qnorm(1-beta/2))/sqrt(t)))

  # library(ldbounds)
  # probcuta = bounds(t,alpha=alpha)$exit.pr
  # probcutb = bounds(t,alpha=alpha)$exit.pr

  stop.eff=list()
  stop.fut=list()
  cutpoints=matrix(0,length(reviews),3)

  for(i in 1:length(reviews)){
    stop.eff[[i]]=which(pbinom(0:reviews[i],reviews[i],p0)>(1-probcuta[i]))-1
    stop.fut[[i]]=c(-1,which(pbinom(0:reviews[i],reviews[i],p1)<(probcutb[i]))-2)
    cutpoints[i,]=c(reviews[i],max(stop.fut[[i]]),min(stop.eff[[i]]))
    cutpoints[i,2]=ifelse(cutpoints[i,2]==-1,NA,cutpoints[i,2])
    cutpoints[i,3]=ifelse(cutpoints[i,3]>cutpoints[i,1],NA,cutpoints[i,3])
  }
  cutpoints[i,]=c(reviews[i],min(stop.eff[[i]])-1,min(stop.eff[[i]]))
  cutpoints=data.frame(cutpoints)

  names(cutpoints)=c("n","futility","efficacy")

  table=matrix(,6,length(r)+1)
  table[1,1]="Stop early - Futility"
  table[2,1]="Stop early - Efficacy"
  table[3,1]="Continue to final analysis - Efficacy"
  table[4,1]="Continue to final analysis - Futility"
  table[5,1]="Inconclusive"
  table[6,1]="Expected number of patients recruited"




  results=list()
  for(m in 1:length(r)){
    # variable to save the results for early exit
    # set up outcome matrix with initial state 0 patients.
    outcome=rep(0,max(reviews)+1)
    outcome[1]=1

    # prepare tempoury variables for outcome
    noutcome=rep(0,max(reviews)+1)


    results[[m]]=matrix(0,0,3)
    patients=0
    for(i in 1:length(reviews)){


      # simulate from operating characteristics
      for(j in 1:(patients+1)){
        noutcome[j:(j+reviews[i]-patients)]=noutcome[j:(j+reviews[i]-patients)]+outcome[j]*dbinom(0:(reviews[i]-patients),reviews[i]-patients,r[m])
      }
      if(reviews[i]<max(reviews)){
        results[[m]]=rbind(results[[m]],0)
        results[[m]][dim(results[[m]])[1],]=c(reviews[i],2,sum(noutcome[stop.eff[[i]]+1]))
        results[[m]]=rbind(results[[m]],0)
        results[[m]][dim(results[[m]])[1],]=c(reviews[i],1,sum(noutcome[stop.fut[[i]]+1]))
        noutcome[c(stop.eff[[i]]+1,stop.fut[[i]]+1)]=0
      } else {
        results[[m]]=rbind(results[[m]],0)
        results[[m]][dim(results[[m]])[1],]=c(reviews[i],3,sum(noutcome[stop.eff[[i]]+1]))
        results[[m]]=rbind(results[[m]],0)
        results[[m]][dim(results[[m]])[1],]=c(reviews[i],4,sum(noutcome[stop.fut[[i]]+1]))
        results[[m]]=rbind(results[[m]],0)
        results[[m]][dim(results[[m]])[1],]=c(reviews[i],5,sum(noutcome[-c(stop.eff[[i]]+1,stop.fut[[i]]+1)]))
        noutcome=rep(0,max(reviews)+1)
      }

    outcome=noutcome
    noutcome=rep(0,max(reviews)+1)
    patients=reviews[i]
    }
  table[1,m+1]=sum(subset(results[[m]][,3],results[[m]][,2]==1))*100
  table[2,m+1]=sum(subset(results[[m]][,3],results[[m]][,2]==2))*100
  table[3,m+1]=sum(subset(results[[m]][,3],results[[m]][,2]==3))*100
  table[4,m+1]=sum(subset(results[[m]][,3],results[[m]][,2]==4))*100
  table[5,m+1]=sum(subset(results[[m]][,3],results[[m]][,2]==5))*100
  table[6,m+1]=sum(results[[m]][,1]*results[[m]][,3])


  }

  # Prepare tables for output
  data=data.frame("x"=table[,1])
  st=c("Stopping rules")
  for(i in 2:dim(table)[2]){
    data=data.frame(data,x=round(as.numeric(table[,i]),2))
    st=c(st,paste0("R=",r[i-1]))
  }
  data=setnames(data,st)

  alpha=sum(data[2:3,2])
  power=sum(data[2:3,3])
  exp.p0=data[6,2]
  exp.p1=data[6,3]
  eta=min(round(1-pbeta(p0,prior.a+cutpoints$futility+1,prior.b+cutpoints$n-cutpoints$futility-1),3),na.rm = TRUE)
  zeta=min(round(pbeta(p1,prior.a+cutpoints$efficacy,prior.b+cutpoints$n-cutpoints$efficacy),3),na.rm = TRUE)


  print(cutpoints)
  cat("\n")
  print(data)


  return(trialDesign_binom_one(reviews=cutpoints$n,success=cutpoints$efficacy,failure=cutpoints$futility,eta=eta,zeta=zeta,alpha=alpha,power=power,exp.p0=exp.p0,exp.p1=exp.p1,p0=p0,p1=p1))

}
