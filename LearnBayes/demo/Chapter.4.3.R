###################################################
# Section 4.3 A Multinomial Model
###################################################

library(LearnBayes)

 alpha = c(728, 584, 138)
 theta = rdirichlet(1000, alpha)

 hist(theta[, 1] - theta[, 2], main="")

S=readline(prompt="Type  <Return>   to continue : ")

###########################################

data(election.2008)
attach(election.2008)

prob.Obama=function(j)
{
p=rdirichlet(5000,
  500*c(M.pct[j],O.pct[j],100-M.pct[j]-O.pct[j])/100+1)
mean(p[,2]>p[,1])
}

Obama.win.probs=sapply(1:51,prob.Obama)

sim.election=function()
{
winner=rbinom(51,1,Obama.win.probs)  
sum(EV*winner)         
}

sim.EV=replicate(1000,sim.election())

windows()
hist(sim.EV,min(sim.EV):max(sim.EV),col="blue")
abline(v=365,lwd=3)  # Obama received 365 votes
text(375,30,"Actual \n Obama \n total")
