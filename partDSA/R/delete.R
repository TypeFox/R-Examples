####################################################################
###                                                                #
### Delete is a function which takes all available basis functions #
###  and looks at combining two possible functions                 #
###                                                                #
####################################################################

check.del <- function(bas.fx, real.LIST, y, wt, opts) {
  try.del <- best.bfs.to.combine(bas.fx=bas.fx, y=y, wt=wt, opts)
  if(!is.matrix(bas.fx))
    bas.fx <- as.matrix(bas.fx)
  
  ## make new.del.bf with try.del
  if(!is.null(try.del$bf.1) ) {
    new.del.list <- get.corresponding.list(bfs=try.del, real.LIST=real.LIST)
    ##is it clean?
  } else {
    new.del.list <- NULL  # if only one bf - return all empty items
  }
  
  return(list(new.bf=new.del.list, del.risk=try.del$del.risk, newbas.fx=(bas.fx[,try.del$bf.1]|bas.fx[,try.del$bf.2])))
}

### return number of the two best basis functions to combine
best.bfs.to.combine <- function(bas.fx, y, wt, opts) {
  if(!is.matrix(bas.fx))
    bas.fx <- as.matrix(bas.fx)

  if(ncol(bas.fx) > 1) {
    M <- ncol(bas.fx) # number of basis functions
    all.combo <- matrix(NA, nrow=M-1, ncol=M)
    for(i in 1:(M-1)) {
      for(j in (i+1):M) {
        new.bf <- bas.fx[,i] + bas.fx[,j]
        new.bf <- ifelse(new.bf > 1, 1, new.bf)
        remove.bf <- cbind(bas.fx[,-c(i,j)], new.bf)
        all.combo[i,j]<-risk.fx(bas.fx=data.frame(remove.bf),y=y,wt=wt,opts=opts)
      }
    }
    index.min <- min(all.combo, na.rm=TRUE)
    pick.m <- which.min(all.combo)[1]
    bf.2 <- ceiling(pick.m / (M-1))
    bf.1 <- pick.m - ((M-1) * (bf.2-1))
  } else {
    bf.1 <- bf.2 <- index.min <- NULL
  }

  return(list(bf.1=bf.1, bf.2=bf.2, del.risk=index.min))
}

### make new list with by combining two basis functions
get.corresponding.list <- function(bfs, real.LIST) {
  bf1 <- bfs$bf.1
  bf2 <- bfs$bf.2
  num.or.1 <- length(real.LIST[[bf1]][[1]])
  num.or.2 <- length(real.LIST[[bf2]][[1]])
  new.bf <- real.LIST
  for(P in 1:length(real.LIST[[bf1]])) {
    for(K in 1:num.or.2) {
      new.bf[[bf1]][[P]][[num.or.1 + K]] <- real.LIST[[bf2]][[P]][[K]]
    }
  }
  new.bf[bf2] <- NULL
  new.bf[[bf1]] <- check.for.repeated.OR.statments(new.list=new.bf[[bf1]])
  new.bf[[bf1]] <- check.for.contiguous.OR.statments(new.list=new.bf[[bf1]])
  return(new.bf)
}

check.for.contiguous.OR.statments <- function(new.list) {
  p <- length(new.list)
  k <- length(new.list[[1]])

  ## combine any regions which are listed twice (in the same bf)
  K <- 1
  while(K < k) { # hold one K and compare to subsequent k's
    compare.mat <- matrix(NA, nrow=p, ncol=(k-K))
    for (P in 1:p) {
      k.mat <- matrix(unlist(new.list[[P]]), byrow=T, ncol=2)
      compare.mat[P,] <- ifelse(apply(matrix(apply(k.mat, 2, func1, ind.K=K), ncol=2),
                                      1, sum) == 2, 1, 0)
      ## compare.mat is a p by (k-K)
    }
    keep <- apply(compare.mat, 2, sum) == (p - 1)

    if(sum(as.numeric(keep)) > 0) {
      STOP <- 0
      get.next.col <- 1
      while(!STOP) {
        K.col <- c(1:length(keep))[keep] # list of which columns are different
        if(length(K.col) > 1)
          K.col <- K.col[get.next.col]

        ## which variable are the k's not the same?
        P.row <- c(1:p)[compare.mat[,K.col] == 0]
        first.K <- K
        second.K <- K + K.col
        reg1 <- new.list[[P.row]][[first.K]]
        reg2 <- new.list[[P.row]][[second.K]]

        contig <- NULL
        if(reg1[1] == reg2[2])
          contig <-c(reg2[1], reg1[2])
        if(reg2[1]==reg1[2])
          contig <- c(reg1[1], reg2[2])

        ## if there are contiguous regions - combine and then
        ## start again at outer while loop - with K=1
        if(!is.null(contig)) {
          new.list[[P.row]][[first.K]] <- contig
          for(P in 1:p)
            new.list[[P]][second.K] <- NULL
          k <- length(new.list[[1]])
          K <- 1
          STOP <- 1
        }

        ## if this isn't contiguous and there are no more regions
        ## to look at - go to the next K
        if(is.null(contig) && get.next.col == sum(keep)) {
          K <- K + 1
          STOP <- 1
        }

        ## if this isn't contiguous but other k's might be
        ## with K redo the inner while loop
        if(is.null(contig) && get.next.col < sum(keep)) {
          get.next.col <- get.next.col + 1
        }

      }
    } else {
      K <- K + 1
    }
  }

  return(new.list)
}

func1<-function(mat,ind.K){
  return(mat[ind.K] == mat[(ind.K+1):length(mat)])      
}

check.for.repeated.OR.statments <- function(new.list) {
  p <- length(new.list)  # number of variables
  k <- length(new.list[[1]])  # number of or statements

  ## combine any regions which are listed twice (in the same bf)
  K <- 1
  while(K < k) {
    compare.mat <- matrix(NA, nrow=p, ncol=(k-K))
    for (P in 1:p) {
      k.mat <- matrix(unlist(new.list[[P]]), byrow=TRUE, ncol=2)
      
      compare.mat[P,] <- ifelse(apply(matrix(apply(k.mat, 2, func1, ind.K=K), ncol=2),
                                      1, sum) == 2, 1, 0)
      if(length(which(compare.mat[P,]==1))>0 & k.mat[K,1]==-Inf & k.mat[K,2]==Inf)
        compare.mat[P,which(compare.mat[P,]==1)]<- NA
    }

    keep <- apply(compare.mat, 2, sum,na.rm=TRUE)
    if(length(which(keep>0)) > 0 & length(which(keep==apply(!is.na(compare.mat),2,sum)))>0 ) {
      identify.col<-which(keep>0)
      for(or.col in identify.col){
        identify.repeated.or<-or.col+1
        var.num<-which(compare.mat[,or.col]>0)
        for(var.row in var.num ){
          new.list[[var.row]][identify.repeated.or][[1]] <- c(-Inf,Inf)
          info('found repeat')
        }
        ## k <- length(new.list[[1]])
      }
    } else {
      K <- K + 1
    }
  }
  return(new.list)
}
