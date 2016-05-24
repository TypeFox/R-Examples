
## Posterior predictive checks under the bottleneck model for the italian data

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

## set to the path to the location of the programs "ms" and "sample_stats" on your computer
## if they are in the current directory you don't need to do anything
pathtoms <- "/"#"/your/path/to/ms/"

## number of loci
numloc<-50
## number of years per generation
numgen<-25
## number of simulations for the posterior check
simpost<-1000
## mutation rate
mutrate <- 2.5*10^(-8)
## recombination  rate
c <- rlnorm(numloc*simpost, -18.148, .5802^2)
## number of base pairs
L<-2000

stat.italy.sim <- stat.3pops.sim[100001:150000,] ## select the simulations for the bottleneck model

## perform the abc for inferring the parameters
par.post<-abc(target= stat.voight["italian",], param=par.italy.sim, sumstat=stat.italy.sim, tol=.005, hcorr = TRUE,method="neuralnet", transf = c("logit","logit","logit","logit"),logit.bounds=rbind(c(0,30000),c(1,100),c(2500,10000),c(40000,60000)))

## sampling with replacement in the multivariate posterior distribution
newsamp<-sample(1:(dim(par.post$adj)[1]),size=simpost,replace=T,prob=par.post$weights) 

newsamp<-par.post$adj[newsamp,]

## duration and start of the bottleneck in the right unit
par1<-(newsamp[,4]-newsamp[,3])/(numgen*4*newsamp[,1])
par2<-(newsamp[,4])/(numgen*4*newsamp[,1])


par.post<-data.frame(theta=rep(newsamp[,1]*4*mutrate*L,each=numloc),
                       rho=rep(newsamp[,1],each=numloc)*4*c*(L-1),
                       L=rep(L, numloc*1000),
                       time.rec=rep(par1,each=numloc),
                       recov=rep(1/newsamp[,2],each=numloc),
                       time.bot=rep(par2,each=numloc),
                       bott=rep(1,times=1000*numloc))
                       
write.table(par.post, file="bottpost", quote=F, row.names=F, col.names=F)
print("Running ms for the bottleneck model...")
### ms simulation and comutations of sum stats
system(paste(".", pathtoms, "ms 10", 1000*numloc, "-t tbs -r tbs tbs -eN tbs tbs -eN tbs tbs < bottpost |", ".", pathtoms, "sample_stats > ss_postb.txt", sep=""))

post.bott <- read.table("ss_postb.txt")
###Average the sum stats over the 50 loci
post.bott <- data.frame(pi.m=tapply(post.bott[,2]/2000, factor(rep(1:1000, each=numloc)), mean,na.rm=T),
                       ss.m=tapply(post.bott[,4], factor(rep(1:1000, each=numloc)), mean,na.rm=T),
                       TajD.m=tapply(post.bott[,6], factor(rep(1:1000, each=numloc)), mean,na.rm=T),
                       TajD.v=tapply(post.bott[,6], factor(rep(1:1000, each=numloc)), var,na.rm=T),
                       pi=tapply(post.bott[,2]/2000, factor(rep(1:1000, each=numloc)), var,na.rm=T),
                       ss.v=tapply(post.bott[,4], factor(rep(1:1000, each=numloc)), var,na.rm=T))

post.bott<-post.bott[,c(1,3,4)]

## by running this line you can re-generate the "ppc" data of the package
save(post.bott, file="ppc.rda", compress=T)
