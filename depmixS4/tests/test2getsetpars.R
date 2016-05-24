# 
# test getpars and setpars
# 

require(depmixS4)

data(speed)

trstart=c(0.896,0.104,0.084,0.916)
trstart=c(trstart[1:2],0,0.01,trstart[3:4],0,0.01)
instart=c(0,1)
resp <- c(5.52,0.202,0.472,0.528,6.39,0.24,0.098,0.902)

mod <- depmix(list(rt~1,corr~1),data=speed,family=list(gaussian(),multinomial()),transition=~Pacc,trstart=trstart,instart=instart,respst=resp,nst=2,ntimes=c(168,134,137))

mod1 <- setpars(mod,getpars(mod))

cat("Test getpars and setpars: ", all.equal(getpars(mod),getpars(mod1)), "(getpars and setpars) \n")


