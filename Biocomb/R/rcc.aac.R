
#########
#1)
#########
pareto.front<- function(scores, classes){
  if (length(scores)!=length(classes)){return(print("scores and classes must be the same length"))}
  if (length(classes)!=length(classes[classes==0])+length(classes[classes==1])){return(print("class values must be 0 or 1"))}


  sorted.scores <- sort(scores,index.return=T)
  ones.at.score <- c()
  zeros.at.score <- c()

  scores.no.duplicates <- c(sorted.scores$x[1])
  if (classes[sorted.scores$ix[1]]==0){
    ones.at.score <- c(0)
    zeros.at.score <- c(1)
  }else{
    ones.at.score <- c(1)
    zeros.at.score <- c(0)
  }

  last.score <- sorted.scores$x[1]
  for (i in 2:length(sorted.scores$x)){
    if (last.score<sorted.scores$x[i]){
      scores.no.duplicates <- c(scores.no.duplicates,sorted.scores$x[i])
      if (classes[sorted.scores$ix[i]]==0){
        ones.at.score <- c(ones.at.score,0)
        zeros.at.score <- c(zeros.at.score,1)
      }else{
        ones.at.score <- c(ones.at.score,1)
        zeros.at.score <- c(zeros.at.score,0)
      }
      last.score <- sorted.scores$x[i]}else{
        if (classes[sorted.scores$ix[i]]==0){
          zeros.at.score[length(zeros.at.score)] <- zeros.at.score[length(zeros.at.score)]+1}else{
            ones.at.score[length(ones.at.score)] <- ones.at.score[length(ones.at.score)]+1
          }
      }
  }

  increased.max.score <- scores.no.duplicates[length(scores.no.duplicates)] + 0.1
  scores.no.duplicates <- c(scores.no.duplicates,increased.max.score)
  ones.at.score <- c(ones.at.score,0)
  zeros.at.score <- c(zeros.at.score,0)

  change.at.score.zero <- !(ones.at.score==0)
  change.at.score.zero[length(change.at.score.zero)] <- T

  change.at.score.one.a <- !(ones.at.score>0 & zeros.at.score==0)
  change.at.score.one <- c(T,change.at.score.one.a[1:(length(change.at.score.one.a)-1)])

  change.at.score <- change.at.score.zero & change.at.score.one

  dup.temp<-scores.no.duplicates
  for (i in 2:length(scores.no.duplicates))
  {scores.no.duplicates[i]<-(dup.temp[i-1]+dup.temp[i])/2  }

  tp <- c()
  tn <- c()
  fp <- c()
  fn <- c()

  pareto.scores <- c()

  tp.next <- sum(classes)
  tn.next <- 0
  fp.next <- length(classes) - tp.next
  fn.next <- 0

  for (i in 1:(length(change.at.score)-1)){
    if (change.at.score[i]){
      pareto.scores <- c(pareto.scores,scores.no.duplicates[i])
      tp <- c(tp,tp.next)
      tn <- c(tn,tn.next)
      fp <- c(fp,fp.next)
      fn <- c(fn,fn.next)
    }
    tp.next <- tp.next - ones.at.score[i]
    tn.next <- tn.next + zeros.at.score[i]
    fp.next <- fp.next - zeros.at.score[i]
    fn.next <- fn.next + ones.at.score[i]
  }
  if (change.at.score[length(change.at.score)]){
    pareto.scores <- c(pareto.scores,scores.no.duplicates[length(change.at.score)])
    tp <- c(tp,tp.next)
    tn <- c(tn,tn.next)
    fp <- c(fp,fp.next)
    fn <- c(fn,fn.next)
  }

  threshold <- c(pareto.scores[length(pareto.scores)])
  cost.factor <- c()
  tp.at.threshold <- c(tp[length(tp)])
  tn.at.threshold <- c(tn[length(tn)])
  fp.at.threshold <- c(fp[length(fp)])
  fn.at.threshold <- c(fn[length(fn)])

  if (length(pareto.scores)==1){cost.factor<-Inf}

  ind <- length(pareto.scores)
  while (ind>1){
    cmin <- Inf
    next.ind <- ind
    for (i in 1:(ind-1)){
      cik <- (fp[i]-fp[ind])/(fn[ind]-fn[i])
      if (cik < cmin & cik > 0){
        cmin <- cik
        next.ind <- i
      }
    }
    if (cmin!=Inf){
      ind <- next.ind
      threshold <- c(threshold,pareto.scores[ind])
      cost.factor <- c(cost.factor,cmin)
      tp.at.threshold <- c(tp.at.threshold,tp[ind])
      tn.at.threshold <- c(tn.at.threshold,tn[ind])
      fp.at.threshold <- c(fp.at.threshold,fp[ind])
      fn.at.threshold <- c(fn.at.threshold,fn[ind])
    }
  }
  cost.factor[length(threshold)]=Inf
  return(data.frame(t=threshold,FP=fp.at.threshold,FN=fn.at.threshold,c=cost.factor))
}
#########
#2)
#########
datRCC<-function(dat,attrs, pos.class){
  datRCC<-matrix(ncol=length(attrs)+1,nrow=nrow(dat))

  for(i in 1:length(attrs)){
    datRCC[,(i+1)]<-10^dat[,attrs[i]]
  }

  for(i in 1:nrow(datRCC)){
    if(grepl(pos.class,dat[i,ncol(dat)])){
      datRCC[i,1]<-0
    }else{
      datRCC[i,1]<-1
    }
  }

  colnames(datRCC)<-c("classes",colnames(dat)[attrs])

  return(datRCC)
}

cost.curve<-function(data, attrs.no, pos.Class, AAC=TRUE, n=101, add=FALSE,xlab="log2(c)",ylab="relative costs",
                     main="RCC",lwd=2,col="black",xlim=c(-4,4), ylim=(c(20,120)))
{

  medians=FALSE
  dat.rcc<-datRCC(data, attrs.no, pos.Class)

  scores=dat.rcc[,2]
  classes=dat.rcc[,1]
  minc <- 2^(xlim[1])
  maxc <- 2^(xlim[2])

  at.threshold<-pareto.front(scores,classes)
  cost.naive.switch <- (length(classes)-sum(classes))/sum(classes)

  cost.function <- function(x){
    tc <- exp(x*log(2))
    denom <- length(classes)-sum(classes)
    if (tc < cost.naive.switch){denom <- tc*sum(classes)}
    for (i in 1:length(at.threshold$c)){
      if (tc<=at.threshold$c[i]){return(100*(at.threshold$FP[i]+tc*at.threshold$FN[i])/denom)}
    }  }
  cost.function.v <- Vectorize(cost.function)



if (medians==TRUE){
    dt<-data.frame(scores,classes)
    t1<-median(dt$scores[dt$classes==0])
    t2<-median(dt$scores[dt$classes==1])

    ind<-which.min(abs(at.threshold$t-t1))
    if (t1<at.threshold$t[length(at.threshold$t)-1]){
        maxc<-2*at.threshold$c[length(at.threshold$t)-1]
        } else if (at.threshold$t[ind]<=t1) {
        maxc<-at.threshold$c[ind]-(t1-at.threshold$t[ind])*(at.threshold$c[ind]-at.threshold$c[ind-1])/(at.threshold$t[ind-1]-at.threshold$t[ind])
        } else {
        maxc<-at.threshold$c[ind+1]-(t1-at.threshold$t[ind+1])*(at.threshold$c[ind+1]-at.threshold$c[ind])/(at.threshold$t[ind]-at.threshold$t[ind+1])
    }


    ind<-which.min(abs(at.threshold$t-t2))
    if (t2>at.threshold$t[1]){
        minc<-0.0001
        } else if (at.threshold$t[ind]<=t2) {
        minc<-at.threshold$c[ind]-(t2-at.threshold$t[ind])*(at.threshold$c[ind]-at.threshold$c[ind-1])/(at.threshold$t[ind-1]-at.threshold$t[ind])
        } else {
         minc<-at.threshold$c[ind+1]-(t2-at.threshold$t[ind+1])*(at.threshold$c[ind+1]-at.threshold$c[ind])/(at.threshold$t[ind]-at.threshold$t[ind+1])
    }

}

  c<-seq(log(minc,base=2),log(maxc,base=2),length=n)
  costs<-cost.function.v(c)

  if (add==FALSE){plot(c,costs,type="l",xlim=xlim, ylim=ylim, xlab=xlab,ylab=ylab, main=main,col=col,lwd=lwd)}
  if (add==TRUE){lines(c,costs,col=col,lwd=lwd)}

  if(AAC==TRUE){
    idx = 2:length(c)
    int<- as.double( (c[idx] - c[idx-1]) %*% (costs[idx] + costs[idx-1]))/2
    return(((c[length(c)]-c[1])*100-int)/((c[length(c)]-c[1])*100))
  }
}

aac.value<-function(rcc){
  return(round(1000*rcc)/1000)
}
