trioSim <-
function(n,popfs,hapfs,edists,recomb=0,riskmod,batchsize=1000){
  # n: number of trios to sample
  # popfs: population frequencies
  # hapfs: list with as many haplotype distributions as there are populations
  # edists: list of functions with as many E-distributions as populations
  # recomb: recombination rate between causal and test loci
  # riskmod: function to calculate disease risk
  # batchsize: Repeat following until n affected trios: 
  #              Generate batchsize trios & extract affecteds

  require(gtools)
  tdat<-NULL
  # Continue sampling batches of trios until we have at least n of them.
  # The test for the number of trios is not nrow(tdat), which gives an
  # error when tdat is NULL.
  while(length(tdat[,1])<n) {
    # sample batchsize trios and add to tdat. Continue until we have
    # at least n trios.
    tdat<-rbind(tdat,tsbatch(popfs,hapfs,edists,recomb,riskmod,batchsize))
  }
  # Convert tdat to a data frame and name columns as expected by Jean/JEM's
  # simulation functions. Can also trim tdat to first n trios.
  tdat<-data.frame(tdat[1:n,])
  names(tdat)<-c("parent1","parent2","child","subpop","attr")
  return(tdat)
}

tsbatch <-
function(popfs,hapfs,edists,recomb,riskmod,batchsize) {
  # Sample batchsize trios.

  tdat<-NULL # initialize
  subpops<-(1:length(popfs))
  # First sample sub-population for trios.
  S<-sample(subpops,batchsize,replace=TRUE,prob=popfs)
  # Loop over the sub-populations. For the trios in each sub-population,
  # sample informative parental haplotype pairs (Hp), child haplotype 
  # pairs (H), E, and D. Return the affected trios.
  for(i in subpops) {
    ind<-(S==i)
    nsub<-sum(ind)
    Hpsub<-samHp(hapfs[[i]],nsub) ##
    Hsub<-samHfromHp(Hpsub,nsub,recomb)
    Esub<-samE(edists[[i]],nsub)
    Dsub<-samD(Hsub,Esub,riskmod)

    dind<-(Dsub==1)
    if(any(dind)) {
      # put together information on the trios. NB: sub-population 
      # indexing by Jean/JEM code is 0,1,... so have to subtract 1 off
      # of i to get sub-pop label.
      tsub<-cbind(HptoGp.p(Hpsub[dind,,drop=FALSE]),
                  HtoGp(Hsub[dind,,drop=FALSE]),i-1,Esub[dind])
      tdat<-rbind(tdat,tsub)
     } 
  }
  return(tdat)
}


samHp <-
function(hapf,n) {
  # FOR verifying samHp function
  # Sample informative parental haplotype pairs.
  # Sample first parent conditional on it being heterozygous, then
  # can sample second parent without any restriction.
  #
  # 1. sample first parent cond'l on het
  # Assign numeric labels 1,2,3,4 to haplos N0, N1, R0, R1, resp.
  pairs<-cbind(rep(1:4,each=4),rep(1:4,times=4))
  probs<-hapf[pairs[,1]]*hapf[pairs[,2]]
  joint.prob = rep(NA,(16*16))
  for(i in 0:15){joint.prob[(16*i+1):(16+16*i)]=probs*probs[(i+1)] }
  joint.id = matrix(1:16^2,nrow=16, ncol=16)
  hetpar.index = matrix(NA,nrow=16*16, ncol=2)
  hetpar.index[,1] <- as.vector(mapply(rep, 1:16, 16))
  hetpar.index[,2] <- rep(c(1:16),16)
  
  hetpairs.id <- c(2,4,5,7,10,12,13,15)
  hetpar1.id <- joint.id[hetpairs.id,]
  hetpar2.id <- joint.id[,hetpairs.id]
  joint.hetpar <- unique(c(hetpar1.id,hetpar2.id)) ## rid of doubly counted
  Hp.id <- sample(joint.hetpar,n,replace=TRUE,prob=joint.prob[joint.hetpar])
  Hp <- hetpar.index[Hp.id,]
  Hp1 <- pairs[Hp[,1],]
  Hp2 <- pairs[Hp[,2],]
  return(cbind(Hp1,Hp2))
}

samHfromHp <-
function(Hp,n,recomb) {
  # Mate the two parents to produce the offspring haplotype pair.
  # Recall: Hp is a 4-column matrix. First pair of columns is the 
  # first parent's haplotype pair, second pair is for second parent.

  if(recomb!=0) stop("non-zero recombination not yet implemented")
  h1s<-sample(1:2,n,replace=TRUE) # choose a haplo of first parent
  H1<-Hp[,1]
  H1[h1s==2]<-Hp[h1s==2,2]
  h2s<-sample(3:4,n,replace=TRUE) # choose a haplo of second parent
  H2<-Hp[,3]
  H2[h2s==4]<-Hp[h2s==4,4]
  return(cbind(H1,H2))
}

samE <-
function(edist,n) {
  # Sample E. Just use the random deviate generator in edist.

  E<-edist(n)
  return(E)
}

samD <-
function(H,E,riskmod) {
  # Sample D. Just need to extract G from H and use the risk model.

  G<-HtoG(H)
  D<-riskmod(G,E)
  return(D)
}

## for riskmod() function


HtoGp <-
function(H) {
  #Extract G' from H
  Gp<-rep(0,nrow(H)) # initialize
  ind<-odd(H[,1]+H[,2]) # trick: heterozyous G' when h1+h2 is odd (i.e., one is odd, the other is even)
  Gp[ind]<-1
  ind<-even(H[,1])&even(H[,2]) # trick: 2 copies when h1 & h2 both even
  Gp[ind]<-2
  return(Gp)
}

HtoG <-
function(H) {
  #Extract G from H, a 2-column matrix of haplotype pairs
  G<-rep(1,nrow(H)) #initialize
  ind<-H[,1]<=2 & H[,2]<=2 # h1 & h2 1 or 2 means 0 copies of R
  G[ind]<-0
  ind<-H[,1]>=3 & H[,2]>=3 # h1 & h2 3 or 4 means 2 copies of R
  G[ind]<-2
  return(G)
}

HptoGp.p <-
function(Hp) {
  # Extract G'_p from Hp
  # JS: added 'drop=FALSE' for Hp with nrow(Hp)=1
  Gp.p<-cbind(HtoGp(Hp[, 1:2, drop=FALSE]),HtoGp(Hp[, 3:4, drop=FALSE]))  
  return(Gp.p)
}


