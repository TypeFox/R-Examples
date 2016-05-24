IsomaxT <- function(x, y, niter){

x.res <- x
dat.mat <- y

obs.di <-  IsoGenem(x.res,dat.mat)
dir.obs <- obs.di[[11]]

n.permut <- niter
n.tests <- nrow(dat.mat)
wf.counts.E <-wf.counts.W<-wf.counts.WC<-wf.counts.M <-wf.counts.I<- matrix(0,n.tests,1)
obs.di.E <- obs.di.W <- obs.di.WC <- obs.di.M <- obs.di.I <- c(1: n.tests)
perm.di.E <- perm.di.W <- perm.di.WC <- perm.di.M <- perm.di.I <- c(1: n.tests)


perm.di <- permut.di <- u.di <- c(1:n.tests)


westfall.E<-westfall.1.E<-rep(0,n.tests)
westfall.W<-westfall.1.W<-rep(0,n.tests)
westfall.WC<-westfall.1.WC<-rep(0,n.tests)
westfall.M<-westfall.1.M<-rep(0,n.tests)
westfall.I<-westfall.1.I<-rep(0,n.tests)


obs.di.E[dir.obs=="u"]<-(obs.di[[1]][dir.obs=="u"])
obs.di.E[dir.obs=="d"]<-(obs.di[[6]][dir.obs=="d"])
obs.di.E1 <- sort(abs(obs.di.E))

obs.di.W[dir.obs=="u"]<-(obs.di[[2]][dir.obs=="u"])
obs.di.W[dir.obs=="d"]<-(obs.di[[7]][dir.obs=="d"])
obs.di.W1 <- sort(abs(obs.di.W))

obs.di.WC[dir.obs=="u"]<-(obs.di[[3]][dir.obs=="u"])
obs.di.WC[dir.obs=="d"]<-(obs.di[[8]][dir.obs=="d"])
obs.di.WC1 <- sort(abs(obs.di.WC))


obs.di.M[dir.obs=="u"]<-(obs.di[[4]][dir.obs=="u"])
obs.di.M[dir.obs=="d"]<-(obs.di[[9]][dir.obs=="d"])
obs.di.M1 <- sort(abs(obs.di.M))

obs.di.I[dir.obs=="u"]<-(obs.di[[5]][dir.obs=="u"])
obs.di.I[dir.obs=="d"]<-(obs.di[[10]][dir.obs=="d"])
obs.di.I1 <- sort(abs(obs.di.I))


pb <- tkProgressBar(title = "progress bar", min = 0, max = niter, 
        width = 300)

for(i in 1 : n.permut){


setTkProgressBar(pb, i, title = paste(round(i/niter * 
                100, 0), "% done"))


xiter.index <- t(sapply(1:n.permut,function(i) sample(x.res)))

perm.di <- IsoGenem(xiter.index[i,],dat.mat)
dir.perm <- perm.di[[11]]

perm.di.E[dir.perm=="u"]<-(abs(perm.di[[1]][dir.perm=="u"]))
perm.di.E[dir.perm=="d"]<-(abs(perm.di[[6]][dir.perm=="d"]))
perm.di.E1 <- sort(perm.di.E)

perm.di.W[dir.perm=="u"]<-(abs(perm.di[[2]][dir.perm=="u"]))
perm.di.W[dir.perm=="d"]<-(abs(perm.di[[7]][dir.perm=="d"]))
perm.di.W1 <- sort(perm.di.W)

perm.di.WC[dir.perm=="u"]<-(abs(perm.di[[3]][dir.perm=="u"]))
perm.di.WC[dir.perm=="d"]<-(abs(perm.di[[8]][dir.perm=="d"]))
perm.di.WC1 <- sort(perm.di.WC)


perm.di.M[dir.perm=="u"]<-(abs(perm.di[[4]][dir.perm=="u"]))
perm.di.M[dir.perm=="d"]<-(abs(perm.di[[9]][dir.perm=="d"]))
perm.di.M1 <- sort(perm.di.M)

perm.di.I[dir.perm=="u"]<-(abs(perm.di[[5]][dir.perm=="u"]))
perm.di.I[dir.perm=="d"]<-(abs(perm.di[[10]][dir.perm=="d"]))
perm.di.I1 <- sort(perm.di.I)



##for E2
permut.di<-perm.di.E[order(obs.di.E)]
u.di[1]<-permut.di[1]

for(j in 2:n.tests){
u.di[j]=max(permut.di[c(1:j)])
#u.mat[j,i]<-max(permut.di[c(1:j)])
}

 test2=(u.di>=obs.di.E1)
 wf.counts.E=test2+wf.counts.E



##for Williams
permut.di<- perm.di.W[order(obs.di.W)]
u.di[1]<-permut.di[1]

for(j in 2:n.tests){
u.di[j]=max(permut.di[c(1:j)])
#u.mat[j,i]<-max(permut.di[c(1:j)])
}

 test2=(u.di>=obs.di.W1)
 wf.counts.W=test2+wf.counts.W



##for Marcus'
permut.di<-perm.di.WC[order(obs.di.WC)]
u.di[1]<-permut.di[1]

for(j in 2:n.tests){
u.di[j]=max(permut.di[c(1:j)])
#u.mat[j,i]<-max(permut.di[c(1:j)])
}

 test2=(u.di>=obs.di.WC1)
 wf.counts.WC=test2+wf.counts.WC



##for M
permut.di<-perm.di.M[order(obs.di.M)]
u.di[1]<-permut.di[1]

for(j in 2:n.tests){
u.di[j]=max(permut.di[c(1:j)])
#u.mat[j,i]<-max(permut.di[c(1:j)])
}

 test2=(u.di>=obs.di.M1)
 wf.counts.M=test2+wf.counts.M



##for Modified M
permut.di<-perm.di.I[order(obs.di.I)]
u.di[1]<-permut.di[1]

for(j in 2:n.tests){
u.di[j]=max(permut.di[c(1:j)])
#u.mat[j,i]<-max(permut.di[c(1:j)])
}

 test2=(u.di>=obs.di.I1)
 wf.counts.I=test2+wf.counts.I

#cat(paste(i,".",sep=""))

}


westfall.E=wf.counts.E/niter
westfall.W=wf.counts.W/niter
westfall.WC=wf.counts.WC/niter
westfall.M=wf.counts.M/niter
westfall.I=wf.counts.I/niter

 k= n.tests-1

for(i in k:1 )
  {
 westfall.1.E[n.tests]=westfall.E[n.tests]
{westfall.1.E[i]<-min(max(westfall.E[i:n.tests]),1) }


  westfall.1.W[n.tests]=westfall.W[n.tests]
{westfall.1.W[i]<-min(max(westfall.W[i:n.tests]),1) }



 westfall.1.WC[n.tests]=westfall.W[n.tests]
{westfall.1.WC[i]<-min(max(westfall.WC[i:n.tests]),1) }


  westfall.1.M[n.tests]=westfall.M[n.tests]
{westfall.1.M[i]<-min(max(westfall.M[i:n.tests]),1) }

   westfall.1.I[n.tests]=westfall.I[n.tests]
{westfall.1.I[i]<-min(max(westfall.I[i:n.tests]),1) }


  }

westfall.E2<-westfall.W2<-westfall.WC2<-westfall.M2<-westfall.I2<-NULL


for (i in 1:n.tests){
westfall.E2[order(obs.di.E)[i]]=westfall.1.E[i]
westfall.W2[order(obs.di.E)[i]]=westfall.1.W[i]
westfall.WC2[order(obs.di.E)[i]]=westfall.1.WC[i]
westfall.M2[order(obs.di.E)[i]]=westfall.1.M[i]
westfall.I2[order(obs.di.E)[i]]=westfall.1.I[i]
}


max5 <- data.frame(row.names(y),westfall.E2,westfall.W2,westfall.WC2,westfall.M2,westfall.I2)
colnames(max5) <-c("Probe.ID","E2", "Williams", "Marcus", "M", "ModM")

close(pb)

return(max5)
}