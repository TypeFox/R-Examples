scale.for.ranks=function(type=1, base=2, silent=TRUE, str.disp.mult=1){

#cat("strength differential: win-tie-loss  Goal differential: median, 25%>this, 10%>this\n")
for(spr in c(.1,.5,1,1.5,2,2.5,3,3.5)){
n=100000
higher.team.goals = rpois(n,exp(spr/log(base)))
lower.team.goals = rpois(n,exp(-spr/log(base)))
win=sum(higher.team.goals>lower.team.goals)/n
loss=sum(higher.team.goals<lower.team.goals)/n
tie=sum(higher.team.goals==lower.team.goals)/n
up90.gd=quantile(abs(higher.team.goals-lower.team.goals),prob=.9)
up75.gd=quantile(abs(higher.team.goals-lower.team.goals),prob=.75)
median.gd=median(abs(higher.team.goals-lower.team.goals))
#cat(paste(format(spr,nsmall=1),": ",round(win*100),"-",round(100*tie),"-",round(100*loss),   " GD: median=",median.gd," up25=",up75.gd, " up10=",up90.gd, "\n",sep=""))
}

goals=goals25=c()
spr=seq(-1.5,2,.5)
for(i in spr){
n=1000000
tmp=rpois(n,exp(i/log(base)))
goals = c(goals,median(tmp))
goals25=c(goals25,quantile(tmp,prob=.75))
}
gdat=rbind(goals,goals25)
colnames(gdat)=format(spr,width=4)
rownames(gdat)=c("median","upper 25%")

if(type==1) spr.seq=seq(0,3,.25)
if(type==2) spr.seq=log(c(1:10,20,50))
win=tie=loss=c()
for(spr in spr.seq){
n=100000
#spr is what user sees
#I'm dividing by log(base) to get scale in table
#need to multiply by log(base) to get back to original
#divide by 2 to get approx attack if attack=defense
higher.team.goals = rpois(n,exp(0.5*spr*log(base)))
lower.team.goals = rpois(n,exp(-0.5*spr*log(base)))
win=c(win,sum(higher.team.goals>lower.team.goals)/n)
loss=c(loss,sum(higher.team.goals<lower.team.goals)/n)
tie=c(tie,sum(higher.team.goals==lower.team.goals)/n)
}
if(type==2) spr.seq=exp(spr.seq)

#cat("\nattack minus defense scores and predicted number of goals\n")
#print(gdat,nsmall=0)

dat=rbind(round(100*win),100-round(100*win)-round(100*loss),round(100*loss))
colnames(dat)=format(str.disp.mult*spr.seq,width=4)
rownames(dat)=c("win","tie", "lose")
if(!silent){
  cat("total strength difference and probability of win, tie, or loss\n")
  print(dat,nsmall=0)
  cat("\nExpected number of goals by team A is attack score of team A divided by defense score of team B.\n")
}

invisible(dat)
}


