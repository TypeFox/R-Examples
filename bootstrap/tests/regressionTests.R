library(bootstrap)

bcanon(c(1,2),1000,mean,alpha=c(0.05))$confpoints
bcanon(c(1,2),1000,mean,alpha=c(0.025))$confpoints
bcanon(c(1,2),1000,mean,alpha=c(0.025,0.5,0.75))$confpoints
bcanon(c(1,2),1000,mean,alpha=c(0.01,0.99))$confpoints
