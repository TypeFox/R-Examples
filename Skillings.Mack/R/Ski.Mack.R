#options(repos = c(CRAN="http://cran.r-project.org")) 
Ski.Mack <- function(y, groups=NULL, blocks=NULL,
                    simulate.p.value = FALSE, B = 10000){
  specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
  options(digits=10)
  if(!is.matrix(y)){
    blocks <- factor(blocks)
    k <- nlevels(blocks)
    groups <- factor(groups)
    t <- nlevels(groups)
    d <- matrix(unlist(split(y, blocks)), ncol = k, byrow = FALSE)
    
    if (any(is.na(groups)) || any(is.na(blocks))){
      stop("NA's are not allowed in groups or blocks")}
    if (any(diff(c(length(y), length(groups), length(blocks))))){
      stop("y, groups and blocks must have the same length")}
    for(i in 1:k){
      if(sum(is.na(d[,i])) >= (t-1)){
        stop("Block#",i," has only one observation. Please remove this block")}
    }
  }				
  if(is.matrix(y)){
    d <- as.matrix(y)
    k <- dim(d)[2]
    t <- dim(d)[1]	
  }
  for(i in 1:k){
    if(sum(is.na(d[,i])) >= (t-1)){
      stop("Block#",i," has only one observation. Please remove this block")}
  }
  y.rank.NA <- matrix(NA,nrow=t,ncol=k)
  y.rank <- matrix(NA,nrow=t,ncol=k)
  for (i in 1:k){
    y.rank.NA[,i]<- rank(d[,i], ties.method = c("average"),na.last="keep")
    y.rank[,i]<- ifelse(is.na(y.rank.NA[,i])== TRUE,(sum(!is.na(y.rank.NA[,i]))+1)/2,y.rank.NA[,i]) 
  }                 
  A = matrix(NA,nrow=1,ncol=t)
  for (j in 1:t){
    Ai = matrix(NA,nrow=1,ncol=k)
    for (i in 1:k){
      Ai[1,i] <- sqrt(12/(sum(!is.na(y.rank.NA[,i]))+1))*(y.rank[j,i]-(sum(!is.na(y.rank.NA[,i]))+1)/2)
    }
    A[1,j]<- sum(Ai)
  }
  transform.matrix <- function(mat) {
    num.of.rows <- nrow(mat)
    num.of.cols <- ncol(mat)
    output <- matrix(rep(0, num.of.rows * num.of.rows), c(num.of.rows, num.of.rows))
    for(i in seq(1,num.of.rows-1)) 
      for (j in seq(i+1,num.of.rows)) {
        output[i,j] <- -sum(mat[i,] * mat[j,])
        output[j,i] <- output[i,j]
      }
    for(i in 1:num.of.rows ){
      output[i,i] <- (-1)*sum(output[,i])
    }
    output
  }
  CovMat  <- transform.matrix(!is.na(y.rank.NA))
  Cov.Inv <- ginv(CovMat)
  rank.cov <- matrix.rank(Cov.Inv)
  T <- A%*%Cov.Inv%*%t(A)
  SM.test <- function(y){	
    if(is.matrix(y)){
      d<- as.matrix(y)
      k <- dim(d)[2]
      t <- dim(d)[1]	
    }
    y.rank.NA <- matrix(NA,nrow=t,ncol=k)
    y.rank <- matrix(NA,nrow=t,ncol=k)
    for (i in 1:k){
      y.rank.NA[,i] <- rank(d[,i], ties.method = c("average"),na.last="keep")
      y.rank[,i]    <- ifelse(is.na(y.rank.NA[,i])== TRUE,(sum(!is.na(y.rank.NA[,i]))+1)/2,y.rank.NA[,i]) 
    }
    
    A = matrix(NA,nrow=1,ncol=t)
    for (j in 1:t){
      Ai = matrix(NA,nrow=1,ncol=k)
      for (i in 1:k){
        Ai[1,i] <- sqrt(12/(sum(!is.na(y.rank.NA[,i]))+1))*(y.rank[j,i]-(sum(!is.na(y.rank.NA[,i]))+1)/2)	
      }
      A[1,j]<- sum(Ai)
    }
    transform.matrix <- function(mat) {
      num.of.rows <- nrow(mat)
      num.of.cols <- ncol(mat)
      output <- matrix(rep(0, num.of.rows * num.of.rows), c(num.of.rows, num.of.rows))
      for(i in seq(1,num.of.rows-1)) 
        for (j in seq(i+1,num.of.rows)) {
          output[i,j] <- -sum(mat[i,] * mat[j,])
          output[j,i] <- output[i,j]
        }
      for(i in 1:num.of.rows ){
        output[i,i] <- (-1)*sum(output[,i])
      }
      output
    }
    CovMat   <- transform.matrix(!is.na(y.rank.NA))
    Cov.Inv  <- ginv(CovMat)
    T.SM     <- A%*%Cov.Inv%*%t(A)
  }
  if(simulate.p.value == TRUE){
    logic.d <- is.na(d)
    num.missing <- matrix(NA,nrow=k,ncol=1)
    for(i in 1:k){
      num.missing[i,1] <- sum(logic.d[,i])  
    }
    remove <- function(x,logic){
      x <- x[logic == FALSE]
      return(x)
    }
    simulated.SM <- matrix(NA,nrow=B,ncol=1)

    
    for (b.num in 1:B){
      y.sim <- matrix(NA,nrow=t,ncol=k)
      
      for(i in 1:k){      
        #rand.int<- as.vector(sample(1:(t-num.missing[i,1]),replace = FALSE))
        if(num.missing[i,1] == 0){
          count <- 0
          rand.int<- sample(y.rank[,i],replace = FALSE)
          for (j in 1:t){       
            count = count+1
            y.sim[j,i]<- rand.int[count]   
          }
        }
        if(num.missing[i,1] > 0){
          count <- 0
          rand.int<-  sample(remove(y.rank[,i],logic.d[,i]) ,replace = FALSE) 
          for (j in 1: t){
            if(logic.d[j,i] == TRUE){y.sim[j,i]<- NA}
            if(logic.d[j,i] == FALSE){
              count = count+1
              y.sim[j,i]<- rand.int[count]
            } 
          }
        }   
      }
      simulated.SM[b.num]	<- SM.test(y.sim)
    }    
    sim.pval<- round(sum(simulated.SM >= as.numeric(T))/B,digits = 7)
  }
  if(simulate.p.value == TRUE){
    cat("\n")
    cat(c("Skillings-Mack Statistic = "),specify_decimal(T,6),c(", p-value = "),specify_decimal(pchisq(T, df = rank.cov , lower.tail = FALSE),6),"\n")
    cat(c("Note: the p-value is based on the chi-squared distribution with d.f. = "),rank.cov,"\n")
    cat(c("Based on B = "),specify_decimal(B,0),c(", Simulated p-value = "),specify_decimal(sim.pval,6),"\n")
    cat("\n")
  }
  if(simulate.p.value == FALSE){
    cat("\n")
    cat(c("Skillings-Mack Statistic = "),specify_decimal(T,6),c(", p-value = "),specify_decimal(pchisq(T, df = rank.cov , lower.tail = FALSE),6),"\n")
    cat(c("Note: the p-value is based on the chi-squared distribution with d.f. = "),rank.cov,"\n")
    cat("\n")
  }
  print(list(Nblocks = k, Ntreatments = t, rawdata = d, rankdata =  y.rank,
             varCovarMatrix = CovMat, adjustedSum = A ))
}