.LKK <-
function (case.x1.0,case.x1.1,ctrl.x1.0,ctrl.x1.1,case.x2.0,case.x2.1,ctrl.x2.0,ctrl.x2.1,p1,p5,p6,OR.SNP1,OR.SNP2,OR.ME) {
  odds.p1=p1/(1-p1)                                       #odds of people without SNP1 without SNP2
  odds.p2=odds.p1*OR.SNP2                                 #odds of people without SNP1 with SNP2
  odds.p3=odds.p1*OR.SNP1                                 #odds of people with SNP1 without SNP2
  odds.p4=odds.p1*OR.SNP1*OR.SNP2*OR.ME                   #odds of people with SNP1 with SNP2
  p1=odds.p1/(1+odds.p1)                                  #p1=the prevalence of people without SNP1 without SNP2 (study-level parameter)
  p2=odds.p2/(1+odds.p2)                                  #p2=the prevalence of people without SNP1 with SNP2
  p3=odds.p3/(1+odds.p3)                                  #p3=the prevalence of people with SNP1 without SNP2
  p4=odds.p4/(1+odds.p4)                                  #p4=the prevalence of people with SNP1 with SNP2
  p5=p5                                                   #p5=the MAF of SNP1 in the whole population (study-level parameter)
  p6=p6                                                   #p6=the MAF of SNP2 in population with major allele of SNP1 (study-level parameter)
  P.y=p1*(1-p5)*(1-p6)+p2*(1-p5)*p6+p3*p5*(1-p6)+p4*p5*p6 #P.y=the expected prevalence of whole population

  SNP1.case=(p4*p6+p3*(1-p6))*p5/P.y
  SNP1.ctrl=((1-p4)*p6+(1-p3)*(1-p6))*p5/(1-P.y)
  SNP2.case=(p4*p5+p2*(1-p5))*p6/P.y
  SNP2.ctrl=((1-p4)*p5+(1-p2)*(1-p5))*p6/(1-P.y)

  LKH1=dbinom(case.x1.1,size=case.x1.0+case.x1.1,SNP1.case,log=TRUE)
  LKH2=dbinom(ctrl.x1.1,size=ctrl.x1.0+ctrl.x1.1,SNP1.ctrl,log=TRUE)
  LKH3=dbinom(case.x2.1,size=case.x2.0+case.x2.1,SNP2.case,log=TRUE)
  LKH4=dbinom(ctrl.x2.1,size=ctrl.x2.0+ctrl.x2.1,SNP2.ctrl,log=TRUE)

  return(sum(c(LKH1,LKH2,LKH3,LKH4)))
}
.RW.1 <-
function (p1,p5,p6,step=0.1,p1.range=3) {
  odds.p1=p1/(1-p1)
  odds.p5=p5/(1-p5)
  odds.p6=p6/(1-p6)
  new.odds=exp(rnorm(3,mean=log(c(odds.p1,odds.p5,odds.p6)),sd=c(step*p1.range,step,step)))
  new.odds[new.odds<exp(-10)]=exp(-10)
  new.odds[new.odds>exp(10)]=exp(10)
  return(new.odds/(1+new.odds))
}
.RW.2 <-
function (OR.SNP1,OR.SNP2,OR.ME,step=0.05) {
  new.OR=exp(rnorm(3,mean=log(c(OR.SNP1,OR.SNP2,OR.ME)),sd=c(step,step,step)))
  return(new.OR)
}
.MCMC.step1 <-
function (case.x1.0,case.x1.1,ctrl.x1.0,ctrl.x1.1,case.x2.0,case.x2.1,ctrl.x2.0,ctrl.x2.1,OR.SNP1,OR.SNP2,OR.ME,iterations.step1=20000,step.random.walk=0.1,p1.times=3,progress.bar=TRUE) {
  N.study=length(case.x1.0)
  Outcome.list=list()
  Outcome=matrix(NA,nrow=N.study,ncol=9)
  colnames(Outcome)=c("p1","p1.cil","p1.ciu","p5","p5.cil","p5.ciu","p6","p6.cil","p6.ciu")
  rownames(Outcome)=1:N.study
  se.matrix=matrix(NA,nrow=N.study,ncol=3)
  colnames(se.matrix)=c("p1","p5","p6")
  rownames(se.matrix)=1:N.study
  options(warn=-1)
  if (progress.bar) {pb <- txtProgressBar(max = N.study*iterations.step1, style=3)}
  for (i in 1:N.study) {
    chain = array(dim = c(iterations.step1+1,3))
    chain[1,] = c(0.01,(ctrl.x1.1[i]+0.5)/(ctrl.x1.0[i]+ctrl.x1.1[i]+1),(ctrl.x2.1[i]+0.5)/(ctrl.x2.0[i]+ctrl.x2.1[i]+1))
    for (k in 1:iterations.step1) {
      chain[1+k,]=.RW.1(chain[k,1],chain[k,2],chain[k,3],step.random.walk)
      pre.LKH=.LKK(case.x1.0[i],case.x1.1[i],ctrl.x1.0[i],ctrl.x1.1[i],case.x2.0[i],case.x2.1[i],ctrl.x2.0[i],ctrl.x2.1[i],chain[k,1],chain[k,2],chain[k,3],OR.SNP1,OR.SNP2,OR.ME)
      post.LKH=.LKK(case.x1.0[i],case.x1.1[i],ctrl.x1.0[i],ctrl.x1.1[i],case.x2.0[i],case.x2.1[i],ctrl.x2.0[i],ctrl.x2.1[i],chain[k+1,1],chain[k+1,2],chain[k+1,3],OR.SNP1,OR.SNP2,OR.ME)
      prop=exp(post.LKH-pre.LKH)
      if (is.na(prop)) {
        chain[1+k,]=chain[k,]
        } else {
        if (runif(1) > prop) {chain[1+k,]=chain[k,]}
      }
      if (progress.bar) {setTxtProgressBar(pb, k+(i-1)*iterations.step1)}
      if (round(k/300)==k/300) {
        acceptance=1-mean(duplicated(chain[(k-299):k,]))
        if (acceptance<0.2) {step.random.walk=step.random.walk*0.9}
        if (acceptance>0.3) {step.random.walk=step.random.walk*1.1}
      }
    }
    Outcome[i,c(1,4,7)]=apply(log(chain[-(1:iterations.step1/2),]/(1-chain[-(1:iterations.step1/2),])),2,mean)
    se.matrix[i,]=sqrt(diag(cov(log(chain[-(1:iterations.step1/2),]/(1-chain[-(1:iterations.step1/2),])))))
    Outcome[i,c(2,5,8)]=Outcome[i,c(1,4,7)]-qnorm(0.975)*se.matrix[i,]
    Outcome[i,c(3,6,9)]=Outcome[i,c(1,4,7)]+qnorm(0.975)*se.matrix[i,]
    Outcome[i,]=exp(Outcome[i,])/(1+exp(Outcome[i,]))
  }
  if (progress.bar) {close(pb)}
  options(warn=0)
  Outcome.list[[1]]=Outcome
  Outcome.list[[2]]=log(Outcome/(1-Outcome))[,c(1,4,7)]
  Outcome.list[[3]]=se.matrix
  names(Outcome.list)=c("p.matrix","b.matrix","se.matrix")
  return(Outcome.list)
}
.MCMC.step2 <-
function (case.x1.0,case.x1.1,ctrl.x1.0,ctrl.x1.1,case.x2.0,case.x2.1,ctrl.x2.0,ctrl.x2.1,step1.list,iterations.step2=200000,step.random.walk=0.05,show.plot=TRUE,progress.bar=TRUE) {
  chain = array(dim = c(iterations.step2+1,3))
  chain[1,] = c(1,1,1)
  b.matrix=step1.list[[2]]
  se.matrix=step1.list[[3]]
  N.study=length(case.x1.0)
  options(warn=-1)
  if (progress.bar) {pb <- txtProgressBar(max = iterations.step2, style=3)}
  for (k in 1:iterations.step2) {
    chain[1+k,]=.RW.2(chain[k,1],chain[k,2],chain[k,3],step=step.random.walk)
    pre.LKH=.LKK(case.x1.0,case.x1.1,ctrl.x1.0,ctrl.x1.1,case.x2.0,case.x2.1,ctrl.x2.0,ctrl.x2.1,step1.list[[1]][,1],step1.list[[1]][,4],step1.list[[1]][,7],chain[k,1],chain[k,2],chain[k,3])
    post.LKH=.LKK(case.x1.0,case.x1.1,ctrl.x1.0,ctrl.x1.1,case.x2.0,case.x2.1,ctrl.x2.0,ctrl.x2.1,step1.list[[1]][,1],step1.list[[1]][,4],step1.list[[1]][,7],chain[1+k,1],chain[1+k,2],chain[1+k,3])
    prop=exp(post.LKH-pre.LKH)
    if (is.na(prop)) {
      chain[1+k,]=chain[k,]
      } else {
      if (runif(1) > prop) {chain[1+k,]=chain[k,]}
    }
    if (round(k/300)==k/300) {
      acceptance=1-mean(duplicated(chain[(k-299):k,]))
      if (acceptance<0.2) {step.random.walk=step.random.walk*0.9}
      if (acceptance>0.3) {step.random.walk=step.random.walk*1.1}
    }
    if (progress.bar) {setTxtProgressBar(pb,k)}
  }
  if (progress.bar) {close(pb)}
  options(warn=0)
  if (show.plot) {
  par(mfrow=c(2,3))
  hist(log(chain[-c(1:(iterations.step2/2)),1]),nclass=30, main="Posterior of log(OR.SNP1)", xlab="Mean value = red line" )
  abline(v = mean(log(chain[-c(1:(iterations.step2/2)),1])), col="red" )
  hist(log(chain[-c(1:(iterations.step2/2)),2]),nclass=30, main="Posterior of log(OR.SNP2)", xlab="Mean value = red line" )
  abline(v = mean(log(chain[-c(1:(iterations.step2/2)),2])), col="red" )
  hist(log(chain[-c(1:(iterations.step2/2)),3]),nclass=30, main="Posterior of log(OR.ME)", xlab="Mean value = red line" )
  abline(v = mean(log(chain[-c(1:(iterations.step2/2)),3])), col="red" )
  plot(log(chain[-c(1:(iterations.step2/2)),1]),type="l",main="Chain value of log(OR.SNP1)", xlab="Mean value = red line" , ylab="log(OR.SNP1)")
  abline(h = mean(log(chain[-c(1:(iterations.step2/2)),1])), col="red" )
  plot(log(chain[-c(1:(iterations.step2/2)),2]),type="l",main="Chain value of log(OR.SNP2)", xlab="Mean value = red line" , ylab="log(OR.SNP2)")
  abline(h = mean(log(chain[-c(1:(iterations.step2/2)),2])), col="red" )
  plot(log(chain[-c(1:(iterations.step2/2)),3]),type="l",main="Chain value of log(OR.ME)", xlab="Mean value = red line" , ylab="log(OR.ME)")
  abline(h = mean(log(chain[-c(1:(iterations.step2/2)),3])), col="red" )
  }
  Outcome=list()
  Outcome[[1]]=apply(log(chain[-c(1:(iterations.step2/2)),]),2,mean)
  Outcome[[2]]=cov(log(chain[-c(1:(iterations.step2/2)),]))
  Outcome[[3]]=exp(Outcome[[1]])
  Outcome[[4]]=.LKK(case.x1.0,case.x1.1,ctrl.x1.0,ctrl.x1.1,case.x2.0,case.x2.1,ctrl.x2.0,ctrl.x2.1,step1.list[[1]][,1],step1.list[[1]][,4],step1.list[[1]][,7],Outcome[[3]][1],Outcome[[3]][2],Outcome[[3]][3])
  Outcome[[5]]=chain
  names(Outcome)=c("b","vcov","OR","LKH","chain")
  return(Outcome)
}
.ROUND <-
function(X,digits=0) {
  return(formatC(X,format="f",digits=digits))
}
.ROUND.p <-
function(X,p.digits=4) {
  n.X=length(X)
  p=NULL
  for (i in 1:n.X) {
    if (p.digits==0) {p[i]="<1"} else {
      p[i]=as.character(round(X[i],p.digits))
      if (p[i]=="0") {
        if (p.digits==1) {p[i]="<0.1"} else {
          p[i]="<0."
          for (l in 1:(p.digits-1)) {
            p[i]=paste0(p[i],0)
          }
          p[i]=paste0(p[i],1)
        }
      } else {
        p[i]=.ROUND(as.numeric(p[i]),p.digits)
      }
    }
  }
  return(p)
}
