
## Code that has been used to simulate the human data set

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

## set to the path to the location of the programs "ms" and "sample_stats" on your computer
## if they are in the current directory you don't need to do anything
pathtoms <- "/"#"/your/path/to/ms/"

## parameter declarations
L <- 2000
numsim <- 60000
numloc <- 50
numgen <- 25

mutrate <- 2.5*10^(-8)

c <- rlnorm(numloc*numsim, -18.148, .5802^2)
Ne <- rep(myN<-runif(numsim, 0, 30000), each=numloc)

print("Generating constant population size model...")
## constant pop size model
par.const <- data.frame(theta=Ne*4*mutrate*L,
                        rho=Ne*4*c*(L-1),
                        rep(L, numloc*numsim))
write.table(par.const, file="const", quote=F, row.names=F, col.names=F)
## run ms
system(paste(".", pathtoms, "ms 10", formatC(numsim*numloc,digit=7), " -t tbs -r tbs tbs < const |", ".", pathtoms, "sample_stats > afr-const.txt", sep=""))

print("Generating population expansion model...")
## expansion model
b <- rep(10^(runif(numsim,1,2)), each=numloc)
Ne.pres <- Ne*b
par.exp <- data.frame(theta=Ne.pres*4*mutrate*L,
                      rho=Ne.pres*4*c*(L-1),
                      L=rep(L, numloc*numsim),
                      time.exp=rep(runif(numsim, 40000, 60000), each=numloc)/(numgen*4*Ne.pres),
                      rate=1/b)
write.table(par.exp, file="exp", quote=F, row.names=F, col.names=F)
## run ms
system(paste(".", pathtoms, "ms 10", formatC(numsim*numloc,digit=7), " -t tbs -r tbs tbs -eN tbs tbs < exp | ", ".", pathtoms, "sample_stats > afr-exp.txt", sep=""))

print("Generating population bottleneck model...")
## bottleneck model
aux<-10^(runif(numsim,1,2))
r <- rep(aux, each=numloc)
duration<-rep(runif(numsim,2500,10000), each=numloc)
start<-rep(runif(numsim,40000,60000), each=numloc)
par.bott <- data.frame(theta=Ne*4*mutrate*L,
                       rho=Ne*4*c*(L-1),
                       L=rep(L, numloc*numsim),
                       time.rec=(start-duration)/(numgen*4*Ne),
                       recov=1/r,
                       time.bot=(start)/(numgen*4*Ne),
                       bott=1)
write.table(par.bott, file="bott", quote=F, row.names=F, col.names=F)
timeyears<-data.frame(time_years<-par.bott$time.bot*numgen*4*Ne)
## run ms
system(paste(".", pathtoms, "ms 10", formatC(numsim*numloc,digit=7), " -t tbs -r tbs tbs -eN tbs tbs -eN tbs tbs < bott | ", ".", pathtoms, "sample_stats > afr-bott.txt", sep=""))

## Compute sums stat
afr.const <- read.table("afr-const.txt")
afr.const <- data.frame(pi.m=tapply(afr.const[,2]/2000, factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                        ss.m=tapply(afr.const[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                        D.m=tapply(afr.const[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                        D.v=tapply(afr.const[,6], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                        pi.v=tapply(afr.const[,2]/2000, factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                        ss.v=tapply(afr.const[,4], factor(rep(1:numsim, each=numloc)), var,na.rm=T))
afr.exp <- read.table("afr-exp.txt")
afr.exp <- data.frame(pi.m=tapply(afr.exp[,2]/2000, factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                      ss.m=tapply(afr.exp[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                      D.m=tapply(afr.exp[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                      D.v=tapply(afr.exp[,6], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                      pi.v=tapply(afr.exp[,2]/2000, factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                      ss.v=tapply(afr.exp[,4], factor(rep(1:numsim, each=numloc)), var,na.rm=T))
afr.bott <- read.table("afr-bott.txt")
afr.bott <- data.frame(pi.m=tapply(afr.bott[,2]/2000, factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                       ss.m=tapply(afr.bott[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                       D.m=tapply(afr.bott[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                       D.v=tapply(afr.bott[,6], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                       pi.v=tapply(afr.bott[,2]/2000, factor(rep(1:numsim, each=numloc)), var,na.rm=T),
                       ss.v=tapply(afr.bott[,4], factor(rep(1:numsim, each=numloc)), var,na.rm=T))
models <- rep(c("const", "exp", "bott"), each=numsim)


## Restrict the simu to the first 50000 non NA simulations
bool.const<-which(apply(afr.const,FUN=function(x){sum(is.na(x))==0},MARGIN=1))
afr.const<-afr.const[bool.const[1:50000],]

## Restrict the simu to the first 50000 non NA simulations
bool.exp<-which(apply(afr.exp,FUN=function(x){sum(is.na(x))==0},MARGIN=1))
afr.exp<-afr.exp[bool.exp[1:50000],]

## Restrict the simu to the first 50000 non NA simulations
bool.bott<-which(apply(afr.bott,FUN=function(x){sum(is.na(x))==0},MARGIN=1))
afr.bott<-afr.bott[bool.bott[1:50000],]
myN<-myN[bool.bott[1:50000]]
myreduc<-aux[bool.bott[1:50000]]

stat.voight <- data.frame(pi.m=c(.11, .085, .079)/100,
                          ss.m=c(11.1, 7.1, 6.9),	
                          D.m=c(-.2, .28, .18),
                          D.v=c(.55, 1.19, 1.08))
voight <- c(T,F,T,T,F,F)


## generate file
myseq<-seq(1,numsim*50,by=50)

mydur<-duration[myseq]
mydur<-mydur[bool.bott[1:50000]]


mystart<-start[myseq]
mystart<-mystart[bool.bott[1:50000]]

par.italy.sim <- data.frame(Ne=myN, a=myreduc,duration=mydur,start=mystart)

afr.const <- (afr.const[,voight])
afr.exp <- (afr.exp[,voight])
afr.bott <- (afr.bott[,voight])

## join data
stat.3pops.sim <- rbind(afr.const,afr.exp,afr.bott)
names(stat.3pops.sim) <- c("pi", "TajD.m", "TajD.v")
stat.voight <- stat.voight[,voight]
names(stat.voight) <- c("pi","TajD.m", "TajD.v")
row.names(stat.voight) <- c("hausa", "italian", "chinese")

models <- rep(unique(models), each=50000)

attributes(stat.3pops.sim)$na.action <- NULL ## o/w cv4abc fails!!!
attributes(par.italy.sim)$na.action <- NULL

## running this line you can re-generate the "human" data of the package.
save(stat.voight, stat.3pops.sim, par.italy.sim, models, file="human.rda", compress=T)

