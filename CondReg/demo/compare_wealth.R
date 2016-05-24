#' @rdname covwrappers
#' @importFrom timeSeries getDataPart
#' @export
covSample <- function(x.mat, ...){
  mu <- colMeans(x.mat)
  Sigma <- cov(x.mat)
  list(mu=mu, Sigma=Sigma)
}

#' @rdname covwrappers
#' @importFrom timeSeries getDataPart
#' @export
covCondreg <- function(x.mat, ...){
  mu <- colMeans(x.mat)
  x.mat <- sweep(x.mat,2,mu)
  Sigma <- select_condreg(x.mat, kgrid(20,30))$S
  colnames(Sigma) <- rownames(Sigma) <- colnames(x.mat)
  list(mu=mu, Sigma=Sigma)
}

#' @rdname covwrappers
#' @importFrom BurStFin var.shrink.eqcor
#' @importFrom timeSeries getDataPart
#' @export
covLW <- function(x.mat, ...){
  mu <- colMeans(x.mat)
  Sigma <- var.shrink.eqcor(x.mat, ...)
  list(mu=mu, Sigma=Sigma)
}

#' @rdname covwrappers
#' @importFrom BurStFin factor.model.stat
#' @importFrom timeSeries getDataPart
#' @export
covFM <- function(x.mat, ...){
  mu <- colMeans(x.mat)
  Sigma <- factor.model.stat(x.mat, ...)
  list(mu=mu, Sigma=Sigma)
}

## library(condreg)
## load(system.file('data/simulationdata.Rdata',package='condreg'))
cat('hello\n')
library(BurStFin)

y <- kgrid(10,30)

p <- ncol(R)

M <- 45            ## estimation horizon
H <- 15            ## holding period
K <- nrow(R)       ## investment horizon
reltc <- 30*0.0001  ## transaction cost

comp.s <- (M>p)    ## can sample estimate be computed?
comp.cr <- TRUE
comp.lw <- TRUE

sweek <- 101
eweek <- K-H
weekseq <- seq(sweek, eweek, by=H)

## keep track of values used for each rebalancing
allmu <- matrix(0, length(weekseq), p)

## weights
wts <- list()
wts$s <- matrix(0, length(weekseq)+1, p) 
wts$cr <- matrix(0, length(weekseq)+1, p) 
wts$lw <- matrix(0, length(weekseq)+1, p) 

## last period earnings
earnings <- matrix(0, length(weekseq)+1, p)

## wealth
wth.leng <- (length(weekseq)+1)*H

wth <- list()
wth$s <- rep(0, wth.leng);  wth$s[H] <- 1
wth$cr <- rep(0, wth.leng);  wth$cr[H] <- 1
wth$lw <- rep(0, wth.leng);  wth$lw[H] <- 1

## t=1 are initial states and indexing starts at t=2
t <- 2
for (begweek in weekseq) {

  cat('computing week', begweek, '\n')

  trindx <- begweek-c(M:1)  ## train data indx
  Rtr <- scale(R[trindx,], center=TRUE, scale=FALSE)

  ## sample
  if (comp.s) {
    ## compute estimators
    Ssample <- covSample(Rtr)$S
    ## compute new portfolio weights
    wts$s[t,] <- pfweights(Ssample)
    ## subtract transaction cost from last period wealth
    tc <- transcost(wts$s[t,], wts$s[t-1,],
                    earnings[t-1,], reltc,
                    wth$s[(t-1)*H])
    ## wth$s[(t-1)*H] <- wth$s[(t-1)*H] - tc
  }

  ## condreg
  if (comp.cr) {
    ## compute estimators
    ## soln <- select_condreg(Rtr,y)
    Scondreg <- covCondreg(Rtr)$S
    ## compute new portfolio weights    
    wts$cr[t,] <- pfweights(Scondreg)
    ## subtract transaction cost from last period wealth    
    tc <- transcost(wts$cr[t,], wts$cr[t-1,],
                    earnings[t-1,], reltc,
                    wth$cr[(t-1)*H])
    ## wth$cr[(t-1)*H] <- wth$cr[(t-1)*H] - tc
  }

  ## ledoit-wolf
  if (comp.lw) {
    ## compute estimators
    Slw <- covLW(Rtr)$S
    ## Slw <- var.shrink.eqcor(Rtr)
    ## compute new portfolio weights    
    wts$lw[t,] <- pfweights(Slw)
    ## subtract transaction cost from last period wealth    
    tc <- transcost(wts$lw[t,], wts$lw[t-1,],
                    earnings[t-1,], reltc,
                    wth$lw[(t-1)*H])
    ## wth$lw[(t-1)*H] <- wth$lw[(t-1)*H] - tc
  }
  
  ## compute this holding period earnings
  ret <- rep(1,p)
  for (i in 0:H){
    ret <- ret * (1+R[begweek+i,])
    ## browser()
    if (comp.s) wth$s[(t-1)*H+i+1] <- wth$s[(t-1)*H]* (ret %*% wts$s[t,])
    if (comp.cr) wth$cr[(t-1)*H+i+1] <- wth$cr[(t-1)*H]* (ret %*% wts$cr[t,])
    if (comp.lw) wth$lw[(t-1)*H+i+1] <- wth$lw[(t-1)*H]* (ret %*% wts$lw[t,])
  }
  earnings[t,] <- ret

  t <- t + 1

}

lent <- length(wth[[1]])
lenr <- length(wth)

rmat <- matrix(unlist(wth),lent,lenr) ## return matrix
rmat <- rmat[-c(1:(H-1)),]
rmat <- as.data.frame(rmat)
colnames(rmat) <- names(wth)
rmat <- cbind(date=as.Date(rownames(R)[c((head(weekseq,1)-1):(tail(weekseq,1)+H))]),
              rmat)

## if (require('ggplot2') & require(reshape2) & require(car)){
  
##   asdf <- melt(rmat,'date')
##   colnames(asdf) <- c('Date','Method','Wealth')
##   asdf$variable=recode(asdf$Method,"'s'='Sample';'cr'='CondReg';'lw'='Ledoit-Wolf';")
##   p <- ggplot(asdf, aes(Date,Wealth,group=Method,color=Method))+
##     geom_line()+ylim(0,7)+theme_bw()
##   p
  
## } else {

  plot(c(1,length(wth$s)),range(wth$s,wth$cr,wth$lw),type='n')
  lines(wth$s,lty=1,col='black')
  lines(wth$cr,lty=1,col='red')
  lines(wth$lw,lty=1,col='green')
  legend('topleft',
         c('Sample','Ledoit-Wolf','Condreg'),
         col=c('black','green','red'), lty=1)

## }


