subpop.sim <-
function(n=list(stage1=32,enrich=NULL,stage2=32),effect=list(early=c(0,0),final=c(0,0)),
           outcome=list(early="N",final="N"),control=list(early=NULL,final=NULL),
           sprev=0.5,nsim=1000,corr=0,seed=12345678,select="thresh",
           weight=NULL,selim=NULL,level=0.025,method="CT-SD",sprev.fixed=TRUE,file=""){

# outcome functions
time.out <- function(effect,n,standard=TRUE,method="exponential"){
 t.stat <- log(effect[1])-log(effect)
 if(method=="exponential"){
  n.event <- n*(1-exp(-effect))
 }
 var.tstat <-  4/(n.event[1]+n.event)
 if(standard==TRUE){
  t.stat <- -t.stat/sqrt(var.tstat)
 }
 var.stat <- rep(1,length(effect))
 return(list(t.stat=t.stat[2:length(t.stat)],var.stat=var.stat))
}
normal.out <- function(effect,n){
 effect <- effect[1]-effect
 t.stat <- effect*sqrt(n/2)
 var.tstat <-  rep(1,length(effect))
 return(list(t.stat=t.stat[2:length(t.stat)],var.stat=var.tstat))
}
binary.out <- function(effect,n,standard=TRUE,method="LOR"){
 n.risk <- rep(n,length(effect))
 n.event <- n*effect
 if(length(n.event)!=length(n.risk)){stop("need to set length n.risk = n.event")}
 if(sum(n.risk>n.event)!=length(n.risk)){stop("need to set n.risk > n.event")}
 lor <- log(n.event/(n.risk-n.event))
 t.stat <- lor-lor[1]
 var.tstat <-  1/n.event[1]+1/(n.risk[1]-n.event[1])+1/n.event+1/(n.risk-n.event)
 if(standard==TRUE){
  t.stat <- t.stat/sqrt(var.tstat)
 }
 var.stat <- 1/n.event+1/(n.risk-n.event)
 return(list(t.stat=t.stat[2:length(t.stat)],var.stat=var.stat))
}

# reporting summary
report.1 <- function(n,sprev,nsim,enrich,ran.seed,rule,selim,thresh,method,t.level,ofile) {
 cat("\n",file=ofile,append=TRUE)
 cat("asd: simulations for adaptive seamless designs: v2.0: 11/11/2013","\n",sep="",file=ofile,append=TRUE)
 cat("\n",file=ofile,append=TRUE)
 cat("sample sizes (per arm): sub-pop stage 1 =",as.integer(sprev*n$stage1),": sub-pop stage 1 =",n$stage1,"\n",sep=" ",file=ofile,append=TRUE)
 cat("sample sizes (per arm): sub-pop stage 2 =",as.integer(sprev*n$stage2),": sub-pop stage 2 =",n$stage2,"\n",sep=" ",file=ofile,append=TRUE)
 if(enrich==FALSE){
 cat("sample sizes (per arm): enrichment: sub-pop stage 2 =",n$enrich,"\n",sep=" ",file=ofile,append=TRUE)
 }
 cat("simulations: n =",nsim,", seed =",ran.seed,"and rule =",rule,"with limits =",round(selim,2),"\n",sep=" ",file=ofile,append=TRUE)
 t.level <- paste(as.character(round(100*t.level,1)),"%",sep="")
 cat("method:",method,"and level =",t.level," (one-sided)","\n",sep=" ",file=ofile,append=TRUE)
}
report.2 <- function(out.lab,lab,eff.c,eff.s,eff.f,ofile) {
 cat(out.lab,lab,"control =",eff.c,": sub-pop =",round(eff.s,2),": full-pop =",round(eff.f,2),"\n",sep=" ",file=ofile,append=TRUE)
}
report.3 <- function(ecorr.lab,fcorr.lab,correl,ofile) {
 cat("correlation: early",ecorr.lab,"and final",fcorr.lab,"=",round(correl,2),"\n",sep=" ",file=ofile,append=TRUE)
}
report.4 <- function(zearly,z1,z2,weight,ofile) {
 cat("\n",file=ofile,append=TRUE)
 cat("simulation of test statistics:","\n",sep=" ",file=ofile,append=TRUE)
 cat("expectation early: sub-pop =",round(zearly[1],2),": full-pop =",
                             round(zearly[2],2),"\n",sep=" ",file=ofile,append=TRUE)
 cat("expectation final stage 1: sub-pop =",round(z1[1],2),": full-pop =",
                             round(z1[2],2),"\n",sep=" ",file=ofile,append=TRUE)
 cat("expectation final stage 2: sub-pop only =",round(z2[1],2),": full-pop only =",
                                                round(z2[2],2),"\n",sep=" ",file=ofile,append=TRUE)
 cat("expectation final stage 2, both groups selected: sub-pop =",round(z2[3],2),
                                              ": full-pop =",round(z2[4],2),"\n",sep=" ",file=ofile,append=TRUE)
 cat("weights: stage 1 =",round(sqrt(weight),2),"and stage 2 =",round(sqrt(1-weight),2),"\n",sep=" ",file=ofile,append=TRUE)
 cat("\n")
}


report.5 <- function(sim.res,n,ofile) {
  cat("hypotheses rejected and group selection options at stage 1 (n):","\n",file=ofile,append=TRUE)
  cat(format(" ",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format("Hs",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format("Hf",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
        "\t",format("Hs+Hf",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
        "\t",format("Hs+f",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
        "\t",format("n",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
        "\t",format("n%",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),"\n",file=ofile,append=TRUE)
  cat(format("sub",digits=1,trim=TRUE,justify="left",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[1,1],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[1,2],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[1,3],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[1,4],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[1,5],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(round(100*sim.res$results[1,5]/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),"\n",file=ofile,append=TRUE)
  cat(format("full",digits=1,trim=TRUE,justify="left",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[2,1],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[2,2],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[2,3],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[2,4],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[2,5],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(round(100*sim.res$results[2,5]/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),"\n",file=ofile,append=TRUE)
 cat(format("both",digits=1,trim=TRUE,justify="left",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[3,1],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[3,2],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[3,3],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[3,4],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sim.res$results[3,5],digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(round(100*sim.res$results[3,5]/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),"\n",file=ofile,append=TRUE)
 cat(format("total",digits=1,trim=TRUE,justify="left",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sum(sim.res$results[1:3,1]),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sum(sim.res$results[1:3,2]),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sum(sim.res$results[1:3,3]),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sum(sim.res$results[1:3,4]),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sum(sim.res$results[1:3,5]),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format("-",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),"\n",file=ofile,append=TRUE)
 cat(format("%",digits=1,trim=TRUE,justify="left",scientific=FALSE,nsmall=0,width=6),
        "\t",format(round(100*sum(sim.res$results[1:3,1])/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
        "\t",format(round(100*sum(sim.res$results[1:3,2])/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
        "\t",format(round(100*sum(sim.res$results[1:3,3])/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
        "\t",format(round(100*sum(sim.res$results[1:3,4])/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
        "\t",format(round(100*sum(sim.res$results[1:3,5])/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
        "\t",format("-",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),"\n",file=ofile,append=TRUE)
 all.reject <- 100*(sum(sim.res$results[1:3,1])+sum(sim.res$results[1:3,2])-sum(sim.res$results[1:3,3]))/n
 cat("reject Hs and/or Hf = ",paste(round(all.reject,3),"%",sep=""),"\n",sep=" ",file=ofile,append=TRUE)
}


# validate inputs
if (length(n$stage1)!=1 | length(n$stage2)!=1) {
        stop("Invalid sample size for stages 1 and 2")
    }
n$stage1 <- as.integer(n$stage1); n$stage2 <- as.integer(n$stage2)
if(is.null(n$enrich)==FALSE){n$enrich <- as.integer(n$enrich)}

if (length(corr)==1) {
 if (abs(corr)>1) {
  stop("Correlation must be value between -1 and 1")
 }
} else {
  stop("Correlation must be value length between -1 and 1")
}
outcome.options <- c("N","B","T")
e.outcome <- as.integer(match(outcome$early,outcome.options,-1))
if (e.outcome < 1) {
   stop("Unknown early outcome: current options N, B or T")
}
if(length(effect$early)!=2){
 stop("Early effect: need vector length=2")
}
if(length(effect$final)!=2){
 stop("Final effect: need vector length=2")
}
if (level >= 1 | level <= 0) {
  stop("level must be between 0 and 1")
}

 sel.options <- c("thresh","futility")
 isel <- as.integer(match(select, sel.options, -1))
 if (isel < 1) {
   stop("Unknown method: current option thresh")
  } # end if
  if(select=="thresh"){
   if(is.null(selim)==TRUE){
    selim[1] <- -6
    selim[2] <- 6
   }
   if(selim[1]>selim[2]){
     stop("Limits for threshold rule: selim[1]<selim[2]")
   } # end if
  } # end if
  selection.rules <- c("threshold","futility")

 meth.options <- c("CT-Simes","CT-Bonferroni","CT-SD","CEF")
    imeth <- as.integer(match(method, meth.options, -1))
    if (imeth < 1) {
        stop("unknown method: current options CT-Simes, CT-Bonferroni, CT-SD or CEF")
    }
 test.method <- c("CT-Simes","CT-Bonferroni","CT-SD","CEF")


# early outcome
if(outcome$early=="N"){
  if(is.null(control$early)){
   early.S <- c(0,effect$early[1])
   early.F <- c(0,effect$early[2])
  } else {
   early.S <- c(control$early,effect$early[1])
   early.F <- c(control$early,effect$early[2])
  }
  early.out <- normal.out
  oearly.lab <- "normal, standardized effect sizes:"
  ecorr.lab <- "standardized effect"
} else if(outcome$early=="B"){
   lbin.test <- sum(effect$early<=0); ubin.test <- sum(effect$early>=1)
   if(lbin.test>0 | ubin.test>0){
    stop("early effect: rate must be values between 0 and 1")
   }
  if(is.null(control$early)==TRUE){
    stop("control effects: rates must be values between 0 and 1")
  }
  if(control$early<0 | control$early>1){
    stop("control effects: rates must be values between 0 and 1") 
  } 
  early.S <- c(control$early,effect$early[1])
  early.F <- c(control$early,effect$early[2])
  early.out <- binary.out
  oearly.lab <- "binary, event rate:"
  ecorr.lab <- "log(odds)"
} else if(outcome$early=="T"){
  lhaz.test <- sum(effect$early<=0)
  if(lhaz.test>0){
    stop("Early effect: hazard must be > 0")
  }
  if(is.null(control$early)){
   early.S <- c(1,effect$early[1])
   early.F <- c(1,effect$early[2])
  } else {
   early.S <- c(control$early,effect$early[1])
   early.F <- c(control$early,effect$early[2])
  }
  early.out <- time.out
  oearly.lab <- "time-to-event, hazard rate:"
  ecorr.lab <- "log(hazard)"
}

# final outcome
f.outcome <- as.integer(match(outcome$final,outcome.options,-1))
if (f.outcome < 1) {
   stop("unknown early outcome: current options N, B or T")
}
if(length(effect$final)!=length(effect$early)){
 stop("final effect: should be vector of same length as early effect")
}

if(outcome$final=="N"){
  if(is.null(control$final)){
   final.S <- c(0,effect$final[1])
   final.F <- c(0,effect$final[2])
  } else {
   final.S <- c(control$final,effect$final[1])
   final.F <- c(control$final,effect$final[2])
  }
  final.out <- normal.out
  ofinal.lab <- "normal, standardized effect sizes:"
  fcorr.lab <- "standardized effect"
} else if(outcome$final=="B"){
   lbin.test <- sum(effect$final<=0); ubin.test <- sum(effect$final>=1)
   if(lbin.test>0 | ubin.test>0){
    stop("Final effect: rate must be value length between 0 and 1")
   }
   if(is.null(control$final)==TRUE){
    stop("control effects: rates must be values between 0 and 1")
   }
   if(control$final<0 | control$final>1){
    stop("control effects: rates must be values between 0 and 1") 
   } 
   final.S <- c(control$final,effect$final[1])
   final.F <- c(control$final,effect$final[2])
   final.out <- binary.out
   ofinal.lab <- "binary, event rate:"
   fcorr.lab <- "log(odds)"
} else if(outcome$final=="T"){
   lhaz.test <- sum(effect$final<=0)
   if(lhaz.test>0){
    stop("Final effect: hazard must be > 0")
   }
   if(is.null(control$final)){
    final.S <- c(1,effect$final[1])
    final.F <- c(1,effect$final[2])
   } else {
    final.S <- c(control$final,effect$final[1])
    final.F <- c(control$final,effect$final[2])
   }
   final.out <- time.out
   ofinal.lab <- "time-to-event, hazard rate:"
   fcorr.lab <- "log(hazard)"
} 


# set-up simulations
if(is.null(n$enrich)==TRUE){n$enrich <- sprev*n$stage2}

zearly.S <- early.out(effect=early.S,n=sprev*n$stage1)$t.stat
zearly.F <- early.out(effect=early.F,n=n$stage1)$t.stat
zearly <- c(zearly.S,zearly.F)

z1.S <- final.out(effect=final.S,n=sprev*n$stage1)$t.stat
z1.F <- final.out(effect=final.F,n=n$stage1)$t.stat
z1 <- c(z1.S,z1.F)

z2.S1 <- final.out(effect=final.S,n=n$enrich)$t.stat
z2.F1 <- final.out(effect=final.F,n=n$stage2)$t.stat
z2.S2 <- final.out(effect=final.S,n=sprev*n$stage2)$t.stat
z2.F2 <- final.out(effect=final.F,n=n$stage2)$t.stat
z2 <- c(z2.S1,z2.F1,z2.S2,z2.F2)

if(is.null(weight)==TRUE){
 weight <- n$stage1/(n$stage1+n$stage2)
} else if(weight<0 | weight>1){
   stop("weight: must be value between 0 and 1")
}

if(is.null(selim)==TRUE){
  selim[1] <- zearly.S
  selim[2] <- zearly.F
}


# run simulations
if(sprev.fixed==TRUE){

 results <- gsubpop.sim(z.early=zearly,z1=z1,z2=z2,sprev=rep(sprev,2),selim=selim,
            corr=corr,nsim=nsim,seed=seed,level=level,
            select=select,method=method,wt=weight)

} else if(sprev.fixed==FALSE){
s.seed <- runif(n=nsim,min=1,max=10000)

for(k in 1:nsim){

sprev.1 <- rbinom(1,n$stage1,sprev)/n$stage1
sprev.2 <- rbinom(1,n$stage2,sprev)/n$stage2

if(is.null(n$enrich)==TRUE){n$enrich <- sprev.2*n$stage2}

 zearly.S <- early.out(effect=early.S,n=sprev.1*n$stage1)$t.stat
 zearly.F <- early.out(effect=early.F,n=n$stage1)$t.stat
 zearly <- c(zearly.S,zearly.F)
 z1.S <- final.out(effect=final.S,n=sprev.1*n$stage1)$t.stat
 z1.F <- final.out(effect=final.F,n=n$stage1)$t.stat
 z1 <- c(z1.S,z1.F)
 z2.S1 <- final.out(effect=final.S,n=n$enrich)$t.stat
 z2.F1 <- final.out(effect=final.F,n=n$stage2)$t.stat
 z2.S2 <- final.out(effect=final.S,n=sprev.2*n$stage2)$t.stat
 z2.F2 <- final.out(effect=final.F,n=n$stage2)$t.stat
 z2 <- c(z2.S1,z2.F1,z2.S2,z2.F2)

 if(is.null(weight)==TRUE){
  weight <- n$stage1/(n$stage1+n$stage2)
 } else if(weight<0 | weight>1){
   stop("weight: must be value between 0 and 1")
 }

 results.k <- gsubpop.sim(z.early=zearly,z1=z1,z2=z2,sprev=c(sprev.1,sprev.2),selim=selim,
            corr=corr,nsim=1,seed=s.seed[k],level=level,
            select=select,method=method,wt=weight)

 if(k==1){results <- results.k} else {
  results$results <- results$results+results.k$results
 }

}

}


# summarise model
report.1(n=n,sprev=sprev,nsim=nsim,enrich=is.null(n$enrich),ran.seed=seed,rule=selection.rules[isel],
              selim=selim,method=method,t.level=level,ofile=file)
report.2(out.lab="early outcome:",lab=oearly.lab,
          eff.c=round(early.S[1],3),eff.s=round(early.S[2],3),eff.f=round(early.F[2],3),ofile=file)
report.2(out.lab="final outcome:",lab=ofinal.lab,
          eff.c=round(final.S[1],3),eff.s=round(final.S[2],3),eff.f=round(final.F[2],3),ofile=file)
report.3(ecorr.lab,fcorr.lab,correl=corr,ofile=file)

report.4(zearly,z1,z2,weight,ofile=file)
report.5(results,n=nsim,ofile=file)

# output
invisible(results)

}
