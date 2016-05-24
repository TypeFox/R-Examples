ES <-
function(u,maxstop){
  if(missing(maxstop))
{
maxstop <- ncol(u)^2 - ncol(u)
}
 EPS <- sqrt(.Machine$double.eps)
EPS <- EPS
 res <- create.tags(u)
col <- rbind(cbind(res$tags[,1],c(1:nrow(res$tags))*2),cbind(res$tags[,2],c(0:(nrow(res$tags)-1))*2+1))
sortcol <- col[ sort.list(col[,2]),][,1]
col2 <- rbind(cbind(res$tags[,1],c(0:(nrow(res$tags)-1))*2+1),cbind(res$tags[,2],c(1:nrow(res$tags))*2))
sortcol2 <- col2[ sort.list(col2[,2]),][,1]


Xscale <- scaledata(u)$X
Yscale <- scaledata(u)$Y

  m <- (ncol(Xscale))^2 - (ncol(Xscale))
  n <- length(Yscale)

  betaols <- beta <- matrix(0, ncol=m , nrow=m/2+1)
  rsd  <- Yscale
  k <- 1
lambda.step <- NULL
lambda.step2 <- NULL
  max.cor <- matrix(block_cross(Xscale,Yscale,sortcol,sortcol2), nrow=2)
  max.cor <- colSums(max.cor*max.cor)
  jj <- which.max(max.cor)
  ii <- as.vector(rbind(2*jj-1,2*jj))
 max.lambda <- max(abs(matrix(block_cross(Xscale,Yscale,sortcol,sortcol2), nrow=2)[,jj]))
 max.lambda2 <- min(abs(matrix(block_cross(Xscale,Yscale,sortcol,sortcol2), nrow=2)[,jj]))
lambda.step[1] <- max.lambda
lambda.step2[1] <- max.lambda2
  finished <- FALSE
  while(!finished){
     bb <- block_solve(Xscale,Yscale,samplesize=nrow(u),ii,sortcol,sortcol2)

    bcur <- beta[k,ii]
    cjs <- matrix(block_cross(Xscale, rsd,sortcol,sortcol2), nrow=2)
ajs <- matrix(block_cross(Xscale, matrix(block_multiple(Xscale, bb-bcur, ii, sortcol, sortcol2)
,ncol=ncol(u)),sortcol,sortcol2),nrow=2)
    Csq <- colSums(cjs*cjs)
    Asq <- colSums(ajs*ajs)
    P   <- colSums(cjs*ajs)

    CCsq <- Csq[jj][1]
    
    CCsqtAsq <- CCsq-Asq
    CCsqtAsq[jj] <- 0
    gam <- (CCsq-P)
    tt <- sqrt(abs(gam^2 - CCsqtAsq*(CCsq-Csq)))
    gam <- cbind(gam-tt, gam+tt)/CCsqtAsq
jjout <- NULL
if(length(jjout) >= 1)
{
print(length(jj))
jjout <- sort(unique(jjout)) }   
    gam[gam<EPS] <- NA
    gam[jj,] <- NA
if(length(jjout) >= 1)
{
    gam[jjout,] <- NA
}
    step <- pmin(gam[,1], gam[,2], na.rm=TRUE)
if (length(sort(unique(c(jj,jjout))))==maxstop)
{
      step <- 1
}
    if(all(is.na(step))){
      jjnew <- integer()
      step <- 1
      finished <- TRUE
    }else{
      jjnew <- which.min(step)
      step <- min(step, na.rm=TRUE)
    }
    k <- k+1
    beta[k, ii] <- beta[k-1, ii] + step*(bb-bcur)
    betaols[k, ii] <- bb
lambda.step[k] <- (1-step)*lambda.step[k-1]
lambda.step2[k] <- (1-step)*lambda.step2[k-1]
    jj <- c(jj,jjnew)
    ii <- as.vector(rbind(2*jj-1,2*jj))

   rsd <- Yscale - matrix(block_multiple(Xscale, as.matrix(beta[k,]), ii=c(1:length(sortcol)), sortcol, sortcol2),ncol=ncol(u))
if (max(abs(block_cross(Xscale, rsd,sortcol,sortcol2))) < EPS)
{
      finished <- TRUE
}
if (length(sort(unique(c(jj,jjout))))==maxstop)
{
      finished <- TRUE
}
   }
  list(beta=beta, enter.seq=jj, enter.seq.ind=ii, betaols=betaols ,c1=lambda.step, c2=lambda.step2,  maxc= max.lambda)
}
