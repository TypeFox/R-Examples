ETMA <-
function (case.x1.0,case.x1.1,ctrl.x1.0,ctrl.x1.1,case.x2.0,case.x2.1,ctrl.x2.0,ctrl.x2.1,data=NULL,sig.level=0.05,max.step.EM=20,iterations.step1=20000,iterations.step2=200000,start.seed=NULL,show.detailed.plot=TRUE,show.final.plot=TRUE,show.p.matrix=FALSE,progress.bar=TRUE) {
  if (is.null(start.seed)==F) {set.seed(start.seed)}
  dat <- list(x1 = substitute(case.x1.0),
                x2 = substitute(case.x1.1),
                x3 = substitute(ctrl.x1.0),
                x4 = substitute(ctrl.x1.1),
                x5 = substitute(case.x2.0),
                x6 = substitute(case.x2.1),
                x7 = substitute(ctrl.x2.0),
                x8 = substitute(ctrl.x2.1))
  x1 <- eval(dat$x1, data, parent.frame())
  x2 <- eval(dat$x2, data, parent.frame())
  x3 <- eval(dat$x3, data, parent.frame())
  x4 <- eval(dat$x4, data, parent.frame())
  x5 <- eval(dat$x5, data, parent.frame())
  x6 <- eval(dat$x6, data, parent.frame())
  x7 <- eval(dat$x7, data, parent.frame())
  x8 <- eval(dat$x8, data, parent.frame())
  if (length(x1)!=length(x2)|length(x2)!=length(x3)|length(x3)!=length(x4)|length(x4)!=length(x5)|length(x5)!=length(x6)|length(x6)!=length(x7)|length(x7)!=length(x8)) {stop("The lengths of each column are different.\n")} else {
    save.study<-(apply(is.na(cbind(x1,x2,x3,x4,x5,x6,x7,x8)),1,sum)==0)
    x1=x1[save.study]
    x2=x2[save.study]
    x3=x3[save.study]
    x4=x4[save.study]
    x5=x5[save.study]
    x6=x6[save.study]
    x7=x7[save.study]
    x8=x8[save.study]
    if (length(x1)<4) {stop("The number of included studies have to more than 3.\n")} else {
      x1[x1==0]<-1
      x2[x2==0]<-1
      x3[x3==0]<-1
      x4[x4==0]<-1
      x5[x5==0]<-1
      x6[x6==0]<-1
      x7[x7==0]<-1
      x8[x8==0]<-1

      Step1.list=list()
      Step2.list.b=list()
      Step2.list.vcov=list()
      Step2.list.OR=list()
      Step2.list.LKH=list()
      for (a in 1:max.step.EM) {
        if (progress.bar) {cat("Start to step 1 in iteration = ",a,".\n")}
        if (a==1) {
          Step1.list=.MCMC.step1(x1,x2,x3,x4,x5,x6,x7,x8,
                                 OR.SNP1=1,OR.SNP2=1,OR.ME=1,iterations.step1=iterations.step1,progress.bar=progress.bar)
        } else {
          Step1.list=.MCMC.step1(x1,x2,x3,x4,x5,x6,x7,x8,
                                 OR.SNP1=Step2.list.OR[[a-1]][1],OR.SNP2=Step2.list.OR[[a-1]][2],OR.ME=Step2.list.OR[[a-1]][3],iterations.step1=iterations.step1,progress.bar=progress.bar)
        }
        if (progress.bar) {cat("Start to step 2 in iteration = ",a,".\n")}
        Outcome.step2=.MCMC.step2(x1,x2,x3,x4,x5,x6,x7,x8,
                                  step1.list=Step1.list,iterations.step2=iterations.step2,show.plot=show.detailed.plot,progress.bar=progress.bar)
        Step2.list.b[[a]]=Outcome.step2$b
        Step2.list.vcov[[a]]=Outcome.step2$vcov
        Step2.list.OR[[a]]=Outcome.step2$OR
        Step2.list.LKH[[a]]=Outcome.step2$LKH
        cat("Log of likelihood = ",round(Step2.list.LKH[[a]])," in iteration = ",a,".\n\n")
        if (a>=2) {if (round(Step2.list.LKH[[a]])<=round(Step2.list.LKH[[a-1]])) {break}}
      }


      if (show.detailed.plot==FALSE&show.final.plot==TRUE) {
        par(mfrow=c(2,3))
        hist(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),1]),nclass=30, main="Posterior of log(OR.SNP1)", xlab="Mean value = red line" )
        abline(v = mean(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),1])), col="red" )
        hist(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),2]),nclass=30, main="Posterior of log(OR.SNP2)", xlab="Mean value = red line" )
        abline(v = mean(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),2])), col="red" )
        hist(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),3]),nclass=30, main="Posterior of log(OR.ME)", xlab="Mean value = red line" )
        abline(v = mean(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),3])), col="red" )
        plot(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),1]),type="l",main="Chain value of log(OR.SNP1)", xlab="Mean value = red line" , ylab="log(OR.SNP1)")
        abline(h = mean(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),1])), col="red" )
        plot(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),2]),type="l",main="Chain value of log(OR.SNP2)", xlab="Mean value = red line" , ylab="log(OR.SNP2)")
        abline(h = mean(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),2])), col="red" )
        plot(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),3]),type="l",main="Chain value of log(OR.ME)", xlab="Mean value = red line" , ylab="log(OR.ME)")
        abline(h = mean(log(Outcome.step2$chain[-c(1:(iterations.step2/2)),3])), col="red" )
      }
      Final.Outcome=list()
      Final.Outcome[[1]]=Step2.list.b[[a]]
      Final.Outcome[[2]]=Step2.list.vcov[[a]]
      Final.Outcome[[3]]=Step2.list.LKH[[a]]
      Final.Outcome[[4]]=sqrt(diag(Final.Outcome[[2]]))
      Final.Outcome[[5]]=length(x1)-3
      Final.Outcome[[6]]=exp(Final.Outcome[[1]])
      Final.Outcome[[7]]=exp(Final.Outcome[[1]]-qt(1-sig.level/2,Final.Outcome[[5]])*Final.Outcome[[4]])
      Final.Outcome[[8]]=exp(Final.Outcome[[1]]+qt(1-sig.level/2,Final.Outcome[[5]])*Final.Outcome[[4]])
      Final.Outcome[[9]]=Final.Outcome[[1]]/Final.Outcome[[4]]
      Final.Outcome[[10]]=2*pt(abs(Final.Outcome[[9]]),df=Final.Outcome[[5]],lower.tail=FALSE)
      Final.Outcome[[11]]=sig.level
      if (show.p.matrix) {
        Final.Outcome[[12]]=Step1.list$p.matrix
        colnames(Final.Outcome[[12]])=c("p(disease|base)","95%CI.low","95%CI.up","p(SNP1=1)","95%CI.low","95%CI.up","p(SNP2=1)","95%CI.low","95%CI.up")
        names(Final.Outcome)=c("b","vcov","LKK","se","df","OR","ci.l","ci.u","t","pval","sig.level","p.matrix")
      } else {
        names(Final.Outcome)=c("b","vcov","LKK","se","df","OR","ci.l","ci.u","t","pval","sig.level")
      }
      class(Final.Outcome)="ggint"
      return(Final.Outcome)

    }
  }
}


