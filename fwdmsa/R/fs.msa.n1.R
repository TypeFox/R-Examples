fs.MSA.n1 <- function(X, 
                      B,
                      cutoff = default.cutoff,
                      initial.subsample.size = default.initial.subsample.size, 
                      minsize = default.minsize, 
                      seed = default.seed,
                      verbose = TRUE){

#determining n1
#needs functions:
#memberX
#memberR
#X.residual
#check.data

X <- check.data(X)                                                                                              # NEW #
N <- dim(X)[1]
J <- dim(X)[2]
m <- max(X)
data <- list(X=X, N=N, J=J, m=m)

default.minsize <- ifelse(N > 500, floor(N/10), floor(N/5))
default.minsize <- ifelse(N <= 250, floor(N/3), default.minsize)
default.minsize <- ifelse(N <  150, 50, default.minsize)
if (N < minsize) stop("Sample size less than Minsize")

member.Xscore <- memberX(X, minsize=minsize)
member.Rscore <- memberR(X, minsize=minsize)

default.seed <- sample(1:10000,B,replace=F)                                                                     # NEW #
if(length(seed)==1 & length(seed)<B) default.seed <- sample(1:10000,B,replace=F)                                # NEW #
default.initial.subsample.size <- min(member.Rscore$n.mono[1,])*J
default.cutoff <- 5

subsample <- matrix(,N,N)
subsample.multi <- list()
for(b in 1:B){
 set.seed(seed[b])                                                                                             # NEW #
 repeat{
  samp <- sample(1:N, initial.subsample.size)
  tmp.r <- integer()
  for(j in 1:J){
   for(l in 1:member.Rscore$n.mono[1,j]){
    tmp.r <- c(tmp.r,sum(sign(match(member.Rscore$group.member[[j]][[l]], samp,nomatch=0))))
   }
  }
  if(any(tmp.r==0)==FALSE) break
 }
 for(i in length(samp):N){
  subsample[,i] <- c(samp,rep(0, N-length(samp)))
  res.fs <- X.residual(X, samp=samp, member=member.Xscore)
  res <- apply(res.fs^2, 1, sum)
  if(length(samp)==N){break}
  min.samp <- integer()
  for(j in 1:J){
   for(l in 1:member.Rscore$n.mono[1,j]){
    temp1 <- member.Rscore$group.member[[j]][[l]][which(res[member.Rscore$group.member[[j]][[l]]]==min(res[member.Rscore$group.member[[j]][[l]]]))]
    min.samp <- c(min.samp, temp1)
   }
  }
  msamp  <- unique(min.samp)

  for(g in 1:max(member.Xscore)){
   if(sum(member.Xscore[msamp]==g)==0){    
    temp2 <- which(member.Xscore==g)[which(res[which(member.Xscore==g)]==min(res[which(member.Xscore==g)]))]
    min.samp <- c(min.samp, temp2)
   }
  }
  msamp <- unique(min.samp)    
  y <- (1:N)[-msamp]
  samp <- c(msamp, y[order(res[-msamp])][1:(length(samp)+1-length(msamp))])
 }

 subsample.multi[[b]] <- subsample
 if(verbose==TRUE)print(b)
}

# NEW # 
# NEW # 
number.unique.subsample <- integer()
id.unique.subsample <- matrix(,N,B)                                                                               
number.major.subsample <- integer()
for(i in initial.subsample.size:N){
 subsample.member <- matrix(rep(0,N),N,B)
 for(b in 1:B){
  subsample.member[which(1:N %in% subsample.multi[[b]][,i]),b] <- 1
 }
 number.unique.subsample[i] <- length(table(apply(subsample.member,2,paste, collapse="")))
 pattern <- apply(subsample.member,2,paste, collapse="")                                                   
 for(j in 1:number.unique.subsample[i]){                                                                      
  id.unique.subsample[i,which(pattern==names(sort(-table(pattern)))[j])] <- j                             
 }                                                                                                             
 number.major.subsample[i] <- sum(id.unique.subsample[i,]==1)      
} 
n1.unique <- min(which(number.unique.subsample < cutoff))                                                        
n1.major <- min(which((B-number.major.subsample+1) < cutoff))   

fs.output.list <- list(data = data,
                       initial.subsample.size = initial.subsample.size, 
                       subsample.multi = subsample.multi, 
                       number.unique.subsample = number.unique.subsample,                                       
                       number.major.subsample = number.major.subsample,                                             # NEW #
                       id.unique.subsample = id.unique.subsample,                                               # NEW #
                       B = B,
                       seed = seed,                                                                              # NEW # 
                       cutoff = cutoff,
                       n1.unique = n1.unique,                                                                   # NEW #
                       n1.major = n1.major)                                                                         # NEW #
class(fs.output.list) <- "fs.n1.class"
return(fs.output.list)
}
