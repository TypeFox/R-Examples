library(nlme)
library(genefilter)


sigdirPttest <- function(o,alpha=.05){
  adj <- guo(o$unadj,alpha) # significant at level .05?
  sigdir <- ((o$dir*2)-1)*adj
  return(-sigdir)
}

sigdirE2test <- function(o,alpha=.05){
  sigdir <- ((o$dir*2)-1)*(o$adj<alpha)
  return(sigdir)
}

monotonicity <- function(o,alpha=.05){
  sigdir <- sigdirPttest(o,alpha)
  ## no significant trend
  none <- intersect(which((rowSums(sigdir > 0) == 0)),which((rowSums(sigdir < 0) == 0)))
  ## monotonous upward trends 
  up <- intersect(which((rowSums(sigdir > 0) > 0)),which(!(rowSums(sigdir < 0) >
0)))
  ## monotonous downward trends
  down <- intersect(which((rowSums(sigdir < 0) > 0)),which(!(rowSums(sigdir > 0) > 0)))
  ## non-monotonous trends
  anti <- intersect(which((rowSums(sigdir < 0) > 0)),which((rowSums(sigdir > 0) > 0)))
  res <- factor(rep(1:4,length=nrow(sigdir)),labels=c('none','up','down','anti'))
  res[none] <- 'none'
  res[up] <- 'up'
  res[down] <- 'down'
  res[anti] <- 'anti'
  res
}

est.lme <- function(y,ia,ib){
  if(!is.matrix(y)){
    stop('y has to be a numeric matrix')
  }
  B <- nrow(y)
  ng <- length(unique(ia))
  M <- diag(ng)[,ia]
  w <- rowSums(M)
  dM <- t(t(y %*% t(M))/w)
  w <- w/sum(w)
  iMu <- misoreg(dM,w)$result
  iMd <- misoreg(-dM,w)$result
  Hs <- sapply(1:nrow(y),function(n){
    y1 <- y[n,]
    ym <- dM[n,]
    iru <- iMu[n,]
    ird <- iMd[n,]
    RSSu <- sum((y1-iru[ia])^2)
    RSSd <- sum((y1+ird[ia])^2)
    if(RSSd<RSSu) {iru <- ird}
    ik <- rank(iru,ties='min')
    l <- sort(unique(ia))
    astar <- ia
    for(i in 1:length(l)){
      astar[astar==l[i]] <- ik[i]
    }
    iastar <- factor(astar)
    #print(length(levels(iastar)))
    if(length(levels(iastar))==1){
      d <- data.frame(Y=y1,AS=iastar,A=ia,B=ib)
      mylme <- try(lme(Y~1,random=~1|B/A,data=d,control=list(minAbsParApVar=pi/100,returnObject=T)),silent=T)
      if(class(mylme)!='lme') {
        mylme <- try(lme(Y~1,random=~1|B/A,data=d,control=list(minAbsParApVar=exp(0)/200,returnObject=T)),silent=T)
      }
      if(class(mylme)=='lme') {
        v <- VarCorr(mylme)
        Hsg <- as.numeric(v[4,2])^2
        Hsb <- as.numeric(v[2,2])^2
        Hse <- mylme$sigma^2
      } else {
        Hsg <- NA
        Hsb <- NA
        Hse <- NA
      }
     } else {
      d <- data.frame(Y=y1,AS=iastar,A=ia,B=ib)
      mylme <- try(lme(Y~AS,random=~1|B/A,data=d,control=list(minAbsParApVar=pi/100,returnObject=T)),silent=T)
      if(class(mylme)!='lme') {
        mylme <- try(lme(Y~AS,random=~1|B/A,data=d,control=list(minAbsParApVar=exp(0)/200,returnObject=T)),silent=T)
      }
      if(class(mylme)=='lme') {
        v <- VarCorr(mylme)
        Hsg <- as.numeric(v[4,2])^2
        Hsb <- as.numeric(v[2,2])^2
        Hse <- mylme$sigma^2
      } else {
        Hsg <- NA
        Hsb <- NA
        Hse <- NA
      }
 #      cat('rest\n')
    }
    return(c(Hsg,Hsb,Hse))
  })
  Hs <- as.matrix(as.data.frame(Hs))
  rownames(Hs) <- NULL
  colnames(Hs) <- NULL
  Hse <- Hs[3,]
  Hsb <- Hs[2,]
  Hsg <- Hs[1,]
  HICC <-  rep(NA,length(Hse))
  Tss <- rowSums((y-rowMeans(y))^2)
  res <- data.frame(sb=Hsb,sg=Hsg,se=Hse,tss=Tss)
  rownames(res) <- names(Tss)
  return(res)
}




pttest <- function(data,g,B,rep=rep(1,length(g))){
  if(!is.matrix(data)){
    stop('data has to be a numeric matrix')
  }
  if(!is.matrix(B[[1]])){
    l <- sort(unique(g))
    c <- lapply(2:length(l),function(i) 1:length(which(g==l[i-1]|g==l[i])))
    mperms <- lapply(2:length(l),function(i) eperms(B,rep[which(g==l[i-1]|g==l[i])]))
  } else {
    mperms <- B
    B <- nrow(mperms[[1]])
  }
  tt <- srt(data,g,c)
  direction <- tt > 0
  tord <- order(tt,decreasing=T)
  genes <- nrow(data)
  raw <- tt*0
  ptest <- function(n){
    p <- lapply(mperms,function(ps) ps[n,])
    RSS <- srt(data,g,p)
    raw <<- raw+(abs(tt)<=abs(RSS))
    apply(na.omit(tt),2,max)
    NULL
  }
  nDis <- sapply(1:B,ptest)
  p <- (raw+1)/(B+1)
  res <- list(dir=direction,unadj=p)
  return(res)
}

srt <- function(data,g,p){
  ### sequential row t stats
  l <- sort(unique(g))
  c <- lapply(2:length(l),function(i) which(g==l[i-1]|g==l[i]))
  t <- sapply(1:length(c),function(idx) rowttests(data[,c[[idx]]][,p[[idx]]],factor(g[c[[idx]]]),tstatOnly=T)$statistic)
  return(t)
}

guo <- function(p,a=.05){
  pp <- ncol(p)*apply(p,1,min)
  fdr <- p.adjust(pp,method='BH')
  n <- sum(na.omit(fdr)<a)
  res <- p< a*n/(ncol(p)*nrow(p))
  return(res)
}





e2test <- function(data,                    # data (matrix)
                   g,                       # group labels (vector)
                   B,                       # number of permutations or permutation matrix
                   rep=rep(1,length(g))     # random factor grouping variable 
                   ){              
  ## if B is not a permutation matrix generate one
  if(!is.matrix(data)){
    stop('data has to be a numeric matrix')
  }
  if(!is.matrix(B)){
    mperms <- eperms(B,rep)
  } else {
    mperms <- B
    B <- nrow(mperms)
  }
  ## reverse orderd groups
  gd <- (max(g)+1)-g
  ## Total sum of squares
  TSS <- rowSums((data-rowMeans(data))^2)
  ## Sample residual sum of square for upwards/downwards trend 
  tRSSu <- e2(data,g,TSS)
  tRSSd <- e2(data,gd,TSS)
  ## Sample direction
  direction <- tRSSu<tRSSd
  tRSS <- pmin(tRSSu,tRSSd)

  ## number of genes
  genes <- nrow(data)
  ## counting vector global variable that keeps track how often the permutation statistic was more extreme than the original one
  count <- rep(0,genes)
  ## global variable to record raw p-values (if needed)
  raw <- tRSS*0
  
    ## calculate statistic for each permutation and save minimal statistic of each run, increase elements of raw if permutation statistic was more extreme than original

  nDis <- apply(mperms,1,function(p){
    RSS <- pmin(e2(data[,p],g,TSS),e2(data[,p],gd,TSS))
    raw <<- raw+(tRSS>=RSS)
    min(RSS)
  })
  ## how often was the original statistic less extreme than the extremest statistics of each run
  count <- sapply(tRSS,function(t) (sum(t >= nDis)+1))
  res <- list(adj=count/(B+1), # maxT FWE adjusted p.value
              unadj=if(!is.null(raw)) {(raw+1)/(B+1)} else {raw}, # raw p.values, NULL in case of parallel processing 
              dir=direction) # most likely direction of each trend
  return(res)

}

# function that returns a permutation matrix that leaves the variance struction intact


eperms <- function(B,random){
  ro <- order(random)
  rro <- order(ro)
  ms <- c(1,table(random[ro]))
  p <- replicate(B,ro[unlist(sapply(1:(length(ms)-1),function(i) sample(sum(ms[1:i]):(sum(ms[1:(i+1)])-1))))[rro]])
  return(t(p))
}

e2 <- function(data,g,TSS){ # atomic function that calculates the e2 statistic for each row of a matrix
  m <- ncol(data)
  ng <- length(unique(g))
  M <- diag(ng)[,g]
  w <- rowSums(M)
  dM <- t(t(data %*% t(M))/w)
  w <- w/sum(w)
  iM <- misoreg(dM,w)$result 
  iMs <- iM[,g]
  RSS <- rowSums((data-iMs)^2)
  return(RSS/TSS)
}

misoreg <- function (data, weights=rep(1/ncol(data),ncol(data))) {
  .Call("misoreg", PACKAGE = "orQA", data, weights)
}

