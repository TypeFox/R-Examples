BoundModwpt<- function(x, wf="la8", n.levels=4, oldtargets)
{
  
  N <- length(x); storage.mode(N) <- "integer"
  J <- n.levels
  if(2^J > N) stop("wavelet transform exceeds sample size in modwt")
  
  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  ht <- dict$hpf/sqrt(2)
  storage.mode(ht) <- "double"
  gt <- dict$lpf/sqrt(2)
  storage.mode(gt) <- "double"
  
  targets=prepareTargets(oldtargets)#encoding targets using gray code
  
  y <- vector("list", sum(2^(1:J)))
  yn <- length(y)
  crystals1 <- rep(1:J, 2^(1:J))
  crystals2 <- unlist(apply(as.matrix(2^(1:J) - 1), 1, seq, from=0))
  names(y) <- paste("w", crystals1, ".", crystals2, sep="")
  
  W <- numeric(N);V<- numeric(N)
  storage.mode(W)  <- "double"; storage.mode(V)<- "double"
  for(j in 1:J) {
    index <- 0
    jj <- min((1:yn)[crystals1 == j])
    for(n in 0:(2^j / 2 - 1)) {
      index <- index + 1
      #should filter parent node j-1,n
      sc=shouldCompute(c(j-1,n),targets)
      if(sum(sc)>0){
        if(j > 1)
          x <- y[[(1:yn)[crystals1 == j-1][index]]]
        if(n %% 2 == 0) {
          z <- .C("pmodwpt", as.double(x), N, as.integer(j),as.integer(2), L, ht, gt,
                  W = W, V = V, PACKAGE="RHRV")[8:9]
          y[[jj + 2*n + 1]] <- z$W
          y[[jj + 2*n]] <- z$V
        }
        else {
          z <- .C("pmodwpt", as.double(x), N, as.integer(j),as.integer(2), L, ht, gt,
                  W = W, V = V, PACKAGE="RHRV")[8:9]
          y[[jj + 2*n]] <- z$W
          y[[jj + 2*n + 1 ]] <- z$V
        }
      }
      
    }
  }
  
  return(y)
}


shouldCompute2intCode<-function(sc){
  if (codeEquals(sc,c(1,0))){  intCode=0}
  if (codeEquals(sc,c(0,1))) {  intCode=1}
  if (codeEquals(sc,c(1,1))){intCode=2}
  return (intCode)
}

prepareTargets <- function(targets){
  newTargets=list()
  numberTargets=length(targets)/2;
  for (n in 1:numberTargets){ 
    newTargets[[n]]=getC(targets[[2*n-1]],targets[[2*n]])
  }
  return(newTargets)
}

shouldCompute <- function(node,targets){
  
  #compute should be c(0,0)(none), c(1,0)(left),c(0,1)(right), c(1,1) (both)
  compute=c(0,0)
  # node c(0,0) must be computed
  if (nodeEquals(node,c(0,0))){
    len=length(targets)
    for (j in 1:len){
      if (targets[[j]][1]==0){
        compute[[1]]=1;
      }else{
        compute[[2]]=1;
      }
      if (codeEquals(compute,c(1,1))) break  
    }
    
  }else{
    len=length(targets)
    nodeCode=getC(node[[1]],node[[2]])
    codeLen=length(nodeCode)
    compute=c(0,0)
    for (j in 1:len){
      if (length(targets[[j]])>codeLen){
        equals=codeEquals(nodeCode,targets[[j]][1:codeLen])
        if (equals){
          ntc=nodeToCompute(targets[[j]][[codeLen+1]])
          compute=compute+ntc
          compute=compute/max(compute)## avoid (2,0) or (0,2)
          if (codeEquals(compute,c(1,1))) break	
          
        }
      }
    }
    return (compute)
  }
  
  
  
  return(compute);
}


nodeToCompute<-function(number){
  if (number==1){##high pass filter
    return (c(0,1))
  }else{## low pass filter
    return (c(1,0))
  }
}

codeEquals<-function(c1,c2){
  if (length(c1)!=length(c2))
  { 
    return (FALSE)
  }else{
    return (prod(c1==c2)==1)
  }	
}

nodeEquals<-function(n1,n2){
  return ((n1[1]==n2[1])&&(n1[2]==n2[2]));
  
}
