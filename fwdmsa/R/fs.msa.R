fs.MSA <- function(X, 
                   initial.subsample = "random", 
                   initial.subsample.size = default.initial.subsample.size, 
                   minsize = default.minsize, 
                   seed = default.seed,
                   n.low = default.n.low,                                             #NEW# 
                   verbose = TRUE){

#running the forward search for MSA
#needs library(mokken)
#needs functions:
#check.data
#memberX
#memberR
#X.residual
#check.monotonicity.fs
#alpha
#observation.obj.func

X <- check.data(X)								#NEW# 
N <- dim(X)[1]
J <- dim(X)[2]
m <- max(X)
data <- list(X=X, N=N, J=J, m=m)
default.n.low <- N/4

  default.minsize <- ifelse(N > 500, floor(N/10), floor(N/5))
  default.minsize <- ifelse(N <= 250, floor(N/3), default.minsize)
  default.minsize <- ifelse(N <  150, 50, default.minsize)
  if (N < minsize) stop("Sample size less than Minsize")

  default.seed <- round(runif(1,1,10000))

member.Xscore <- memberX(X, minsize=minsize)
member.Rscore <- memberR(X, minsize=minsize)

if(is.numeric(initial.subsample)==FALSE){
 if(initial.subsample=="random"){
   set.seed(seed)
   default.initial.subsample.size <- min(member.Rscore$n.mono[1,])*J
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
 } 
} else {samp <- initial.subsample; initial.subsample.size <- length(initial.subsample)}

if (is.vector(samp)==FALSE) stop("Start sample not defined corretly")

statistic.monitor <- matrix(,J+2,N)
subsample <- matrix(,N,N)
residual.samp <- array(,c(N,J,N))
output.monotonicity.fs <- list()
output.scale.fs <- matrix(,J,N)

for(i in length(samp):N){
    subsample[,i] <- c(samp,rep(0, N-length(samp)))
    res.fs <- X.residual(X, samp=samp, member=member.Xscore)
    res <- apply(res.fs^2, 1, sum)
    residual.samp[,,i] <- res.fs

    output.monotonicity.fs[i] <- list(check.monotonicity.fs(X[samp,], minvi=0.3, minsize=minsize, hi.score=member.Rscore$hi.mono, lo.score=member.Rscore$lo.mono))
    output.scale.fs[,i] <- search.normal(X[samp,], verbose=FALSE)

    statistic.monitor[,i] <- c(alpha(X[samp,]), coefH(X[samp,])$H, coefH(X[samp,])$Hi)
    if(verbose==TRUE)print(i)
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

objective.function <- observation.obj.func(residual.samp, subsample, initial.subsample.size)
order.objective.function <- apply(-objective.function$observation.obj,2,order)


obj.incl <- matrix(,N-1,2)
objective.incl <- matrix(,N,N)
objective.excl <- matrix(,N,N)
for(i in initial.subsample.size:(N-1)){
objective.incl[i,] <- c(objective.function$observation.obj[subsample[,i],i], rep(NA,N-length(objective.function$observation.obj[subsample[,i],i])))
objective.excl[i,] <- c(objective.function$observation.obj[-subsample[,i],i], rep(NA,N-length(objective.function$observation.obj[-subsample[,i],i])))
obj.incl[i,] <- c(boxplot(objective.incl[i,], plot=F)$stats[4]- boxplot(objective.incl[i,],plot=F)$stats[2], boxplot(objective.incl[i,], plot=F)$stats[4])
}
tukey.upper.fence <- cbind(obj.incl[,2]+1.5*obj.incl[,1], obj.incl[,2]+3*obj.incl[,1])
n2.temp <- which((tukey.upper.fence[,1]-objective.function$min.excl)<0)                            # NEW #
n2 <- min(n2.temp[which(n2.temp > n.low)])                                                                    # NEW #
suspect <-  (1:N)[-subsample[,n2]]      							# NEW 


fs.output.list <- list(data=data,
                    initial.subsample.size=initial.subsample.size, 
                    n2=n2,
                    suspect=suspect,
                    subsample=subsample, 
                    member.Xscore=member.Xscore, 
                    member.Rscore=member.Rscore, 
                    residual=residual.samp, 
                    objective=objective.function, 
                    order.objective.function =order.objective.function,
                    tukey.upper.fence=tukey.upper.fence,
                    monotonicity=output.monotonicity.fs, 
                    scale=output.scale.fs, 
                    statistics=statistic.monitor)
class(fs.output.list) <- "fs.class"
return(fs.output.list)
}
