##############################################################################
## Calculate the MMA K1.
K1 <- function(genomat) Kstat(genomat, 1)
##############################################################################
## Calculate the MMA K1+K2.
K12 <- function(genomat) Kstat(genomat, 2)
##############################################################################
## Calculate the MMA K1+K2.
Kstat <- function(genomat, type = 1)
{
  ## calculate the similarity matrix
  similarity <- function(mat){
    nrow <- dim(mat)[1]
    ncol <- dim(mat)[2]
    disMat <- matrix(0,nrow,nrow)
    ones <- matrix(1,ncol,nrow)
    for(i in 1:nrow){
      tmpMat <- abs( mat- t(mat[i,]*ones) )
      disVec <- apply(tmpMat,1,sum)
      disMat[i,] <- disVec
    }  
    return(2*ncol-disMat)
  } ## similarity

  simMat <- similarity(mma.level(genomat))
  result <- sum(simMat) - sum(diag(simMat))
  if(type > 1)
    result <- result + sum(simMat^2) - sum(diag(simMat^2))
  N <- nrow(genomat)
  result / (N * (N - 1))
}
##############################################################################
## standardized genetic dissimilarity--score
eff1 <- function(n, nmark, s1)
{
  sqrt(n) * (2 * nmark - s1) * (n - 1) / (n * nmark)
}

##############################################################################
## function: mma.level
## set 3 levels to 0,1,2
mma.level <- function(mat){
  if(any(is.na(mat)))
    stop("mma routines cannot handle missing values: fill in genotypes")
  as.matrix(mat - min(mat))
}

##############################################################################
## main function: mma-sequentially minimize 1st and 2nd moment
## doesn't handle duplicate rows
## doesn't handle missing values
mma <- function(genof, p, sequent = FALSE, exact = FALSE,
                dismat = FALSE)
{
  nr <- nrow(genof)
  if (p < 1 | p > nr)
    stop(paste("wrong number of rows:", p))
  if (p == 1) {
    warning("no selection when p = 1")
    return(1)
  }
  if (p == nr) {
    warning("all lines selected")
    return(seq(1:nr))
  }
  
  ## function: swap
  ## try to optimize the subset
  swap <- function(disMat,chosen){
    n <- dim(disMat)[1]
    moreSwap <- TRUE

    ## loop until nothing can be swapped
    while(moreSwap) {
      moreSwap <- FALSE

      ## scan rows already chosen
      for(i in 1:n){
        ##cat("i=")
        ##print(i)
        if(chosen[i]==1){
          tmpChosen <- chosen
          tmpChosen[i] <- 0
          tmpMat <- tmpChosen * disMat
          tmpVec <- tmpMat[,i]
          ##score1 <- sum(tmpVec) + sum(tmpVec*tmpVec)
          score1 <- sum(tmpVec)
          
          ## for each already chosen row, find a row not chosen yet
          for(j in 1:n){
            if(chosen[j]==0){
              tmpVec <- tmpMat[,j]
              ##score2 <- sum(tmpVec) + sum(tmpVec*tmpVec)
              score2 <- sum(tmpVec)
              if(score2>score1){
                score1 <- score2
                toSwap <- j
                ##cat("swap")
                ##print(j)
                moreSwap <- TRUE
              }
            }
          }

          if(moreSwap){
            chosen[i] <- 0
            chosen[toSwap] <- 1
            break
          }
        }
      }
      ## cat("out of for loop\n")
    }
    return(chosen)
  } ## swap

  ## function: choose
  ## choose a subset of p
  choose <- function(disMat,p){
    n <- dim(disMat)[1]
    chosen <- rep(0,n)

    ## step1: choose the pair with maximum distance
    maxLoc <- which.max(disMat)
    c1 <- floor(maxLoc / n)
    c2 <- maxLoc - c1*n
    if(c2==0) c2 <- n
    if(c2>0) c1 <- c1 + 1
    chosen[c1] <- 1
    chosen[c2] <- 1
    
    ## step2: add one more row at a time
    if(p==2) {
      return(chosen)
    }

    minVal <- min(disMat)
    for(i in 3:p) { 
      cTmp <- -1
      ## minScore <- n * (maxVal + maxVal^2) + 1
      maxScore <- n * (minVal) 

      ## choose the next row 
      for(j in 1:n) {
        if(chosen[j]==0) {
          tmpVec <- (chosen * disMat)[,j]
          ##score <- sum(tmpVec) + sum(tmpVec*tmpVec)
          score <- sum(tmpVec)

          if(score>maxScore) {
            cTmp <- j
            maxScore <- score
          }
        }
      }

      chosen[cTmp] <- 1
    }
    
    return(chosen)
  } ## choose

  ## function: distance
  ## calculate the distance matrix
  distance <- function(mat, exact = FALSE){
    nrow <- dim(mat)[1]
    ncol <- dim(mat)[2]
    disMat <- matrix(0, nrow, nrow)
    for(i in 1:nrow){
      tmpMat <- abs(mat - matrix(mat[i,], nrow, ncol, byrow = TRUE))
      if(exact)
        tmpMat[tmpMat > 1] <- 1
      disMat[i,] <- apply(tmpMat, 1, sum)
    }  
    disMat
  } ## distance

  ##function op
  op<-function(cList,disMat,cMat,num){
    flag<-1
    update<-0
    while(flag==1){
      flag<-0
      uList<-c(1:ncol(disMat))
      uList<-uList[-cList]
      n<-dim(disMat)[1]

      for( i in 1:num ){
        score1<-sum(disMat[cList[i],cList])-disMat[cList[i],cList[i]] 
        for ( j in 1:(n-num)){
          score2<-sum(disMat[uList[j],cList])-disMat[uList[j],cList[i]]
          if (score1 <= score2) {
            t<-cList
            t[i]<-uList[j]
            if (score1< score2){
              flag<-1 
              ##cat("yes")
              cList<-t
              cMat<-cList
              update<-1
              break
            }
            cMat<-rbind(cMat,t)
          }
        }
        if (flag==1) break
      }
    }
    list(cMat=cMat,update=update)
  }

  ##function op2
  op2<-function(Mat,disMat,num){
    if (is.vector(Mat))
      return (Mat)
    
    tmpMat<-Mat
    workMat<-Mat[-1,]
    tMat<-Mat
    for ( i in 1: dim(Mat)[1]){
      tmpMat[i,]<-sort(tmpMat[i,])
    }
    if(is.vector(workMat)) workMat<-matrix(workMat,1,)
    ##print(dim(workMat)[1])
    while(dim(workMat)[1]>0){
      zeroMat2<-vector(mode="numeric", length=0)
      count<-0
      for ( i in 1: (dim(workMat)[1])){
        zeroMat<-vector(mode="numeric", length=0)
        k<-op(workMat[i,],disMat,zeroMat,num)
        ##cat("in")
        pk<-k$cMat
        update<-k$update
        ##return(pk)
        if(is.vector(pk)) pk<-matrix(pk,1,)
        if (dim(pk)[2]>0){
          if(update==1) {
            tmpMat<-pk
            for ( l in 1: dim(pk)[1]){
              tmpMat[l,]<-sort(tmpMat[l,])
            }
            zeroMat2<-pk
            print(tmpMat)
            count<-count+1
            break
          }
          remove<-rep(1,length=dim(pk)[1])
          for(j in 1: (dim(pk)[1]) ){
            pk[j,]<-sort(pk[j,])       
            test<-abs(t(tmpMat)-pk[j,])
            test1<-apply(test,2,sum)
            if(any(test1==0)) remove[j]<-0
          }
          if (sum(remove)>0){
            count<-count+1
            ##print(sum(remove))
            zeroMat2<-rbind(zeroMat2,pk[remove==1,])
            tmpMat<-rbind(tmpMat,pk[remove==1,])
          }
          ##if (sum(remove)==0){
          ##  zeroMat2<-vector(mode="numeric", length=0)
          ##  dim(zeroMat2)[1]<-0
          ##}
        }
      }
      if (count==0) break
      workMat<-zeroMat2
      if(is.vector(workMat)) workMat<-matrix(workMat,1,)
      ##print(dim(workMat)[1])
      
    }
    tmpMat
  }

  ##function moment2
  moment2<-function(Mat,disMat,num,oriMat,col){
    ## Odd logic to op and op2...
    if(length(Mat) == 1) {
      if(Mat==0)
        return(oriMat$cMat)
    }
    else {
      if(all(Mat == oriMat$cMat))
        return(Mat)
    }

    simMat<- 2*col-disMat
    s2<-rep(0,dim(Mat)[1])
    for(i in 1: dim(Mat)[1]){
      s2[i]<-sum(simMat[Mat[i,],Mat[i,]]^2)
      s2[i]<-(s2[i]-sum(diag(simMat[Mat[i,],Mat[i,]])^2))/2
    } 
    print(s2)
    loc<-which.min(s2)
    print(loc)
    Mat[loc,]
  }
  
  ## calculate the distance matrix
  disMat <- distance(mma.level(genof), exact)

  ## first, choose a subset of p
  chosen <- choose(disMat,p)
  ##print(chosen)
  
  ## second, try to optimize the subset
  chosen <- swap(disMat,chosen)

  ## output the result
  cList <- seq(chosen)[chosen == 1]

  if(sequent) {
    tmp<-op(cList,disMat,cList,p)
    tmp2<-op2(tmp$cMat,disMat,p)
    
    ##cat("Rows chosen:\n")
    ret <- list(cList = cList,
                op = tmp,
                op2 = tmp2,
                moment2 = moment2(tmp2, disMat, p, tmp, dim(genof)[2]))
  }
  else
    ret <- list(cList = cList)
  if(dismat)
    ret$dismat <- disMat
  ret
} ## mma

