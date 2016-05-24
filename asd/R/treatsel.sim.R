treatsel.sim <-
function(n=list(stage1=32,stage2=32),effect=list(early=c(0,0,0),final=c(0,0,0)),
           outcome=list(early="N",final="N"),nsim=1000,corr=0,seed=12345678,select=0,epsilon=1,
           weight=NULL,thresh=1,level=0.025,ptest=c(1),method="invnorm",fu=FALSE,file=""){

# outcome functions
time.out <- function(effect,n,standard=TRUE,method="exponential"){
 t.stat <- log(effect)-log(effect[1])
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
 effect <- effect-effect[1]
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
 t.stat <- lor[1]-lor
 var.tstat <-  1/n.event[1]+1/(n.risk[1]-n.event[1])+1/n.event+1/(n.risk-n.event)
 if(standard==TRUE){
  t.stat <- t.stat/sqrt(var.tstat)
 }
 var.stat <- 1/n.event+1/(n.risk-n.event)
 return(list(t.stat=t.stat[2:length(t.stat)],var.stat=var.stat))
}

# reporting summary
report.1 <- function(ss1,ss2,n,ran.seed,rule,epsilon,thresh,method,t.level,ofile) {
 cat("\n",file=ofile,append=TRUE)
 cat("asd: simulations for adaptive seamless designs: v2.0: 11/11/2013","\n",sep="",file=ofile,append=TRUE)
 cat("\n",file=ofile,append=TRUE)
 cat("sample sizes (per arm): stage 1 =",ss1,": stage 2 =",ss2,"\n",sep=" ",file=ofile,append=TRUE)
 cat("simulations: n =",n,"and seed =",ran.seed,"\n",sep=" ",file=ofile,append=TRUE)
 cat("selection rule:",rule,sep=" ",file=ofile,append=TRUE)
 if(rule=="epsilon rule (select within epsilon of maximum)"){
   cat(" : epsilon =",round(epsilon,3),"\n",file=ofile,append=TRUE)
  } else if(rule=="threshold rule (select greater than or equal to threshold)"){
   cat(" : threshold =",round(thresh,3),"\n",file=ofile,append=TRUE)
  } else {
   cat("\n",file=ofile,append=TRUE)
 }
t.level <- paste(as.character(round(100*t.level,1)),"%",sep="")
cat("method:",method,"and level =",t.level," (one-sided)","\n",sep=" ",file=ofile,append=TRUE)
}
report.2 <- function(out.lab,lab,eff.c,eff.t,ofile) {
 cat(out.lab,lab,"control =",eff.c,": treatment(s) =",eff.t,"\n",sep=" ",file=ofile,append=TRUE)
}
report.3 <- function(ecorr.lab,fcorr.lab,correl,ofile) {
 cat("correlation: early",ecorr.lab,"and final",fcorr.lab,"=",correl,"\n",sep=" ",file=ofile,append=TRUE)
}
report.4 <- function(zearly,z1,z2,weight,ofile) {
 cat("\n",file=ofile,append=TRUE)
 cat("simulation of test statistics:","\n",sep=" ",file=ofile,append=TRUE)
 cat("expectation early =",round(zearly,1),"\n",sep=" ",file=ofile,append=TRUE)
 cat("expectation final stage 1 =",round(z1,1),"and stage 2 =",round(z2,1),"\n",sep=" ",file=ofile,append=TRUE)
 cat("weights: stage 1 =",round(sqrt(weight),2),"and stage 2 =",round(sqrt(1-weight),2),"\n",sep=" ",file=ofile,append=TRUE)
 cat("\n")
}
report.5 <- function(sim.res,n,rej.pow,ofile) {
  cat("number of treatments selected at stage 1:","\n",file=ofile,append=TRUE)
  cat(format(" ",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format("n",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format("%",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),"\n",file=ofile,append=TRUE)
  for (i in 1:length(sim.res$count.total)){
    cat(format(round(i,0),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
     "\t",format(round(sim.res$count.total[i],0),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
     "\t",format(round(100*sim.res$count.total[i]/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
                       "\n",file=ofile,append=TRUE)
    }
  cat(format("Total",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(sum(sim.res$count.total),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format(round(100*sum(sim.res$count.total)/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),"\n",file=ofile,append=TRUE)
 cat("\n",file=ofile,append=TRUE)
 cat("treatment selection at stage 1:","\n",file=ofile,append=TRUE)
  cat(format(" ",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format("n",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format("%",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),"\n",file=ofile,append=TRUE)
  for (i in 1:length(sim.res$count.total)){
    cat(format(round(i,0),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
     "\t",format(round(sim.res$select.total[i],0),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
     "\t",format(round(100*sim.res$select.total[i]/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
              "\n",file=ofile,append=TRUE)
    }
 cat("\n",file=ofile,append=TRUE)
 cat("hypothesis rejection at study endpoint:","\n",file=ofile,append=TRUE)
  cat(format(" ",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format("n",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
        "\t",format("%",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),"\n",file=ofile,append=TRUE)
  for (i in 1:length(sim.res$count.total)){
    cat(format(paste("H",round(i,0),sep=""),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
     "\t",format(round(sim.res$reject.total[i],0),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=6),
     "\t",format(round(100*sim.res$reject.total[i]/n,5),digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=2,width=6),
        "\n",file=ofile,append=TRUE)
    }
 cat("\n",file=ofile,append=TRUE)
 for (i in 1:length(rej.pow)){
  if(i==1){all.hyp <- paste("H",rej.pow[1],sep="")} else if(i>1){
   all.hyp <- append(all.hyp,paste("H",rej.pow[i],sep=""))
  }
 }
 all.hyp <- paste(all.hyp,collapse=" and/or ")
 cat("reject",all.hyp,"=",sim.res$sim.reject,": ",paste(round(100*sim.res$sim.reject/n,3),"%",sep=""),
                  "\n",sep=" ",file=ofile,append=TRUE)
}


# validate inputs
if (length(n$stage1)!=1 | length(n$stage2)!=1) {
        stop("Invalid sample size for stages 1 and 2")
    }
n$stage1 <- as.integer(n$stage1); n$stage2 <- as.integer(n$stage2)
if (length(corr)==1) {
 if (abs(corr)>1) {
  stop("Correlation must be value between -1 and 1")
 }
} else {
  stop("Correlation must be value length equal 1")
}
outcome.options <- c("N","B","T")
e.outcome <- as.integer(match(outcome$early,outcome.options,-1))
if (e.outcome < 1) {
   stop("Unknown early outcome: current options N, B or T")
}
if(length(effect$early)<2){
 stop("Early effect: need vector length > 1")
}
if (level >= 1 | level <= 0) {
  stop("level must be between 0 and 1")
}
ptest.check <- sum(is.element(ptest,seq(1:(length(effect$early)-1))))
if (ptest.check < length(ptest)) {
    stop("invalid ptest vector (see help)")
}
opt.select <- 0:6
e.select <- as.integer(match(select,opt.select,-1))
if (e.select < 1) {
   stop("selection rule: integer betwwen 0 and 6 (see help)")
}
selection.rules <- c("select all treatments","select one","select two",
           "select three","epsilon rule (select within epsilon of maximum)",
           "randomly select a single treatment",
           "threshold rule (select greater than or equal to threshold)")
meth.options <- c("invnorm", "fisher")
    imeth <- as.integer(match(method, meth.options, -1))
    if (imeth < 1) {
        stop("unknown method: current options invnorm or fisher")
    }
comb.method <- c("inverse normal combination test",
           "fisher combination test")


# early outcome
if(outcome$early=="N"){
 if(effect$early[1]!=0){
   stop("Early effect: set value for control treatment to 0")
 }
 early.out <- normal.out
 oearly.lab <- "normal, standardized effect sizes:"
 ecorr.lab <- "standardized effect"

} else if(outcome$early=="B"){
   lbin.test <- sum(effect$early<=0); ubin.test <- sum(effect$early>=1)
   if(lbin.test>0 | ubin.test>0){
    stop("early effect: rate must be values between 0 and 1")
   }
 early.out <- binary.out
 oearly.lab <- "binary, event rate:"
 ecorr.lab <- "log(odds)"

} else if(outcome$early=="T"){
 if(effect$early[1]!=1){
   stop("Early effect: set hazard for control treatment to 1")
 }
  lhaz.test <- sum(effect$early<=0)
  if(lhaz.test>0){
    stop("Early effect: hazard must be > 0")
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
 if(effect$final[1]!=0){
   stop("Final effect: set value for control treatment to 0")
 }
 final.out <- normal.out
 ofinal.lab <- "normal, standardized effect sizes:"
 fcorr.lab <- "standardized effect"

} else if(outcome$final=="B"){
   lbin.test <- sum(effect$final<=0); ubin.test <- sum(effect$final>=1)
   if(lbin.test>0 | ubin.test>0){
    stop("Final effect: rate must be value length between 0 and 1")
   }
 final.out <- binary.out
 ofinal.lab <- "binary, event rate:"
 fcorr.lab <- "log(odds)"

} else if(outcome$final=="T"){
 if(effect$final[1]!=1){
   stop("Final effect: set hazard for control treatment to 1")
 }
  lhaz.test <- sum(effect$final<=0)
  if(lhaz.test>0){
    stop("Final effect: hazard must be > 0")
  }
 final.out <- time.out
 ofinal.lab <- "time-to-event, hazard rate:"
 fcorr.lab <- "log(hazard)"
} 

# set-up for gtreatsel.sim
zearly <- early.out(effect=effect$early,n=n$stage1)$t.stat
vearly <- early.out(effect=effect$early,n=n$stage1)$var.stat
z1 <- final.out(effect=effect$final,n=n$stage1)$t.stat
z2 <- final.out(effect=effect$final,n=n$stage2)$t.stat
v1 <- final.out(effect=effect$final,n=n$stage1)$var.stat
v2 <- final.out(effect=effect$final,n=n$stage2)$var.stat
if(is.null(weight)==TRUE){
 weight <- n$stage1/(n$stage1+n$stage2)
} else if(weight<0 | weight>1){
   stop("weight: must be value between 0 and 1")
}
results <- gtreatsel.sim(z1=z1,z2=z2,zearly=zearly,v1=v1,v2=v2,
         vearly=vearly,corr=corr,weight=weight,epsilon=epsilon,thresh=thresh,
         nsim=nsim,seed=seed,select=select,level=level,ptest=ptest,
         fu=fu,method=method)


# summarise model
report.1(ss1=as.integer(n$stage1),ss2=as.integer(n$stage2),n=nsim,ran.seed=seed,
               rule=selection.rules[e.select],epsilon=epsilon,thresh=thresh,
               method=comb.method[imeth],t.level=level,ofile=file)
report.2(out.lab="early outcome:",lab=oearly.lab,
          eff.c=round(effect$early[1],3),eff.t=round(effect$early[2:length(effect$early)],3),ofile=file)

report.2(out.lab="final outcome:",lab=ofinal.lab,
          eff.c=round(effect$final[1],3),eff.t=round(effect$final[2:length(effect$final)],3),ofile=file)
report.3(ecorr.lab,fcorr.lab,correl=corr,ofile=file)
report.4(zearly,z1,z2,weight,ofile=file)
report.5(results,n=nsim,rej.pow=ptest,ofile=file)

# output
invisible(results)

}
