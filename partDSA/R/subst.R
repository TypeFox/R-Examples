####################################################################
###                                                                #
### Subst is a function which takes all available basis functions  #
###  and looks permutations of all  two possible functions         #
###                                                                #
####################################################################

subst <- function(bas.fx, y, wt, dat, minsplit, minbuck, real.LIST, opts, x.temp,
                  is.num, control) {
  ## first call add - to find best split for each basis function
  try.add <- addon(bas.fx=bas.fx, y=y, wt=wt, dat=dat,
                   minsplit=minsplit, minbuck=minbuck, real.LIST=real.LIST, opts=opts,
                   x.temp=x.temp, is.num=is.num, control=control)

  ## situation 1 - some variables cannot be split
  use.var <- NULL
  for(M in 1:ncol(bas.fx))
    use.var <- c(use.var, try.add[[M]][[4]])
  cant.split <- sum(use.var) >= (ncol(bas.fx) - 1)

  # if there are at least 2 variables to try to do a substution with
  if(!cant.split) {
    ## only grab basis functions which we can use
    bfs.to.split <- c(c(1:length(use.var))[!use.var])
    try.subs <- best.sub <- list()
    best.sub$sub.risk <- Inf
    for(good.to.split in bfs.to.split) {  # come back with best of each one
      other.bfs <- sum(bfs.to.split > good.to.split)
      if(other.bfs > 0) {
        try.subs[[good.to.split]] <- list()
        good.ind <- c(1:length(bfs.to.split))[bfs.to.split == good.to.split]
        for(next.bf in c(bfs.to.split[(good.ind+1):length(bfs.to.split)])) {
          new.regions <- cbind(try.add[[good.to.split]][[2]],
                               try.add[[good.to.split]][[3]],
                               try.add[[next.bf]][[2]],
                               try.add[[next.bf]][[3]])
          try.subs[[good.to.split]][[next.bf]] <-
            bfs.to.combine(bas.fx=new.regions, y=y, wt=wt, opts=opts,
                           var1=good.to.split, var2=next.bf)
          if(try.subs[[good.to.split]][[next.bf]]$sub.risk < best.sub$sub.risk)
            best.sub <- try.subs[[good.to.split]][[next.bf]]
        }
      }
    }
    ## clean up best of subs
    new.sub.list <- get.corresponding.sub.list(bsub=best.sub,
                                               splits=try.add, real.LIST=real.LIST)
  } else {
    ## there are not enough variables to try to do
    ## a subst - make sure to send something back
    ## to indicate that we can't do a split
    best.sub <- new.sub.list <- list()
    best.sub$sub.risk <- NULL
  }
 
  return(list(best.of.subs=new.sub.list, sub.risk=best.sub$sub.risk,
              try.add=try.add,best.sub=best.sub))
} #end subst()

### for substitution only - picks the best of the 6 poss combinations
### return number of the two best basis functions to combine
bfs.to.combine <- function(bas.fx, y, wt, opts, var1, var2) {
  all.poss.subs <- vector("list", 6)
  a <- bas.fx[,1]
  b <- bas.fx[,2]
  c <- bas.fx[,3]
  d <- bas.fx[,4]

  ## a+c b+d
  all.poss.subs[[1]] <- list()
  all.poss.subs[[1]][[2]] <- ifelse(a+c > 1, 1, a+c)
  all.poss.subs[[1]][[3]] <- ifelse(b+d > 1, 1, b+d)
  all.poss.subs[[1]][[1]] <-
        risk.fx(bas.fx=data.frame(cbind(all.poss.subs[[1]][[2]],
                         all.poss.subs[[1]][[3]])),
                         y=y, wt=wt, opts=opts)
  
  ## a+d b+c
  all.poss.subs[[2]] <- list()
  all.poss.subs[[2]][[2]] <- ifelse(a+d > 1, 1, a+d)
  all.poss.subs[[2]][[3]] <- ifelse(b+c > 1, 1, b+c)
  all.poss.subs[[2]][[1]] <-
        risk.fx(bas.fx=data.frame(cbind(all.poss.subs[[2]][[2]],
                         all.poss.subs[[2]][[3]])),
                         y=y, wt=wt, opts=opts)
   
  ## a+b+c d
  all.poss.subs[[3]] <- list()
  all.poss.subs[[3]][[2]] <- ifelse(a+b+c > 1, 1, a+b+c)
  all.poss.subs[[3]][[3]] <- d
  all.poss.subs[[3]][[1]] <-
        risk.fx(bas.fx=data.frame(cbind(all.poss.subs[[3]][[2]],
                         all.poss.subs[[3]][[3]])),
                         y=y, wt=wt, opts=opts)
  
  ## a+b+d c
  all.poss.subs[[4]] <- list()
  all.poss.subs[[4]][[2]] <- ifelse(a+b+d > 1, 1, a+b+d)
  all.poss.subs[[4]][[3]] <- c
  all.poss.subs[[4]][[1]] <-
        risk.fx(bas.fx=data.frame(cbind(all.poss.subs[[4]][[2]],
                         all.poss.subs[[4]][[3]])),
                         y=y, wt=wt, opts=opts)
   
  ## a+c+d b
  all.poss.subs[[5]] <- list()
  all.poss.subs[[5]][[2]] <- ifelse(a+c+d > 1, 1, a+c+d)
  all.poss.subs[[5]][[3]] <- b
  all.poss.subs[[5]][[1]] <-
        risk.fx(bas.fx=data.frame(cbind(all.poss.subs[[5]][[2]],
                         all.poss.subs[[5]][[3]])),
                         y=y, wt=wt, opts=opts)
  
  ## b+c+d a
  all.poss.subs[[6]] <- list()
  all.poss.subs[[6]][[2]] <- ifelse(b+c+d > 1, 1, b+c+d)
  all.poss.subs[[6]][[3]] <- a
  all.poss.subs[[6]][[1]] <-
        risk.fx(bas.fx=data.frame(cbind(all.poss.subs[[6]][[2]],
                         all.poss.subs[[6]][[3]])),
                         y=y, wt=wt, opts=opts)
   
  all.combo <- NULL
  for(ind in 1:6)
    all.combo[ind] <- all.poss.subs[[ind]][[1]]

  index.min <- min(all.combo, na.rm=TRUE)
  pick.sub <- c(1:6)[all.combo == index.min][1]
  bf.2 <- all.poss.subs[[pick.sub]][[3]]
  bf.1 <- all.poss.subs[[pick.sub]][[2]]

  return(list(bf.1.binary=bf.1, bf.2.binary=bf.2, sub.risk=index.min,
              which.sub=pick.sub, bf1=var1, bf2=var2))
}

### make new list with by combining two basis functions
get.corresponding.sub.list <- function(bsub, splits, real.LIST) {
  a <- bf1 <- bsub$bf1  # what is the first bf?
  c <- bf2 <- bsub$bf2  # what is the second bf?

  ## first add on splits
  first.split <- add.bas.fx(BFs=real.LIST, new.add=splits[[bf1]], m=bf1)
  second.split <- add.bas.fx(BFs=real.LIST, new.add=splits[[bf2]], m=bf2)
  both.split.LIST <- add.bas.fx(BFs=first.split, new.add=splits[[bf2]], m=bf2)
  b <- length(real.LIST) + 1
  d <- length(real.LIST) + 2

  ## based on which split is best - now combine (form unions) via deletion
  num.or.a <- length(both.split.LIST[[a]][[1]])
  num.or.b <- length(both.split.LIST[[b]][[1]])
  num.or.c <- length(both.split.LIST[[c]][[1]])
  num.or.d <- length(both.split.LIST[[d]][[1]])
  new.bf <- both.split.LIST

  if(bsub$which.sub == 1) { # a+c b+d
    for(P in 1:length(both.split.LIST[[a]])) {
      for(K in 1:num.or.c) {
        new.bf[[a]][[P]][[num.or.a+K]] <- both.split.LIST[[c]][[P]][[K]]
      }
    }

    for(P in 1:length(both.split.LIST[[b]])) {
      for(K in 1:num.or.d) {
        new.bf[[b]][[P]][[num.or.b+K]] <- both.split.LIST[[d]][[P]][[K]]
      }
    }
    new.bf[[a]] <- check.for.repeated.OR.statments(new.list=new.bf[[a]])
    new.bf[[a]] <- check.for.contiguous.OR.statments(new.list=new.bf[[a]])
    new.bf[[b]] <- check.for.repeated.OR.statments(new.list=new.bf[[b]])
    new.bf[[b]] <- check.for.contiguous.OR.statments(new.list=new.bf[[b]])
    new.bf[c] <- new.bf[d] <- NULL
  }

  if(bsub$which.sub == 2) { # a+d b+c
    for(P in 1:length(both.split.LIST[[a]])) {
      for(K in 1:num.or.d) {
        new.bf[[a]][[P]][[num.or.a+K]] <- both.split.LIST[[d]][[P]][[K]]
      }
    }

    for(P in 1:length(both.split.LIST[[b]])) {
      for(K in 1:num.or.c) {
        new.bf[[b]][[P]][[num.or.b+K]] <- both.split.LIST[[c]][[P]][[K]]
      }
    }

    new.bf[[a]] <- check.for.repeated.OR.statments(new.list=new.bf[[a]])
    new.bf[[a]] <- check.for.contiguous.OR.statments(new.list=new.bf[[a]])
    new.bf[[b]] <- check.for.repeated.OR.statments(new.list=new.bf[[b]])
    new.bf[[b]] <- check.for.contiguous.OR.statments(new.list=new.bf[[b]])
    new.bf[c] <- new.bf[d] <- NULL
  }

  num.or.ab <- length(second.split[[bf1]][[1]])
  num.or.cd <- length(first.split[[bf2]][[1]])

  ## b is the number of bf in real.LIST + 1
  num.or.a <- length(first.split[[a]][[1]])
  num.or.b <- length(first.split[[b]][[1]])
  num.or.c <- length(second.split[[c]][[1]])
  num.or.d <- length(second.split[[b]][[1]])

  if(bsub$which.sub == 3) { # a+b+c d
    new.bf <- second.split
    for(P in 1:length(second.split[[bf1]])) {
      for(K in 1:num.or.c) {
        new.bf[[bf1]][[P]][[num.or.ab+K]] <- second.split[[c]][[P]][[K]]
      }
    }
    new.bf[[bf1]] <- check.for.repeated.OR.statments(new.list=new.bf[[bf1]])
    new.bf[[bf1]] <- check.for.contiguous.OR.statments(new.list=new.bf[[bf1]])
    new.bf[c] <- NULL
  }

  if(bsub$which.sub == 4) { # a+b+d c
    d <- b
    new.bf <- second.split
    for(P in 1:length(second.split[[bf1]])) {
      for(K in 1:num.or.d) {
        new.bf[[bf1]][[P]][[num.or.ab+K]] <- second.split[[d]][[P]][[K]]
      }
    }
    new.bf[[bf1]] <- check.for.repeated.OR.statments(new.list=new.bf[[bf1]])
    new.bf[[bf1]] <- check.for.contiguous.OR.statments(new.list=new.bf[[bf1]])
    new.bf[d] <- NULL
  }

  if(bsub$which.sub == 5) { # a+c+d b
    new.bf <- first.split
    for(P in 1:length(first.split[[bf2]])) {
      for(K in 1:num.or.a) {
        new.bf[[bf2]][[P]][[num.or.cd+K]] <- first.split[[a]][[P]][[K]]
      }
    }
    new.bf[[bf2]] <- check.for.repeated.OR.statments(new.list=new.bf[[bf2]])
    new.bf[[bf2]] <- check.for.contiguous.OR.statments(new.list=new.bf[[bf2]])
    new.bf[a] <- NULL
  }

  if(bsub$which.sub == 6) { # b+c+d a
    new.bf <- first.split
    for(P in 1:length(first.split[[bf2]])) {
      for(K in 1:num.or.b) {
        new.bf[[bf2]][[P]][[num.or.cd+K]] <- first.split[[b]][[P]][[K]]
      }
    }
    new.bf[[bf2]] <- check.for.repeated.OR.statments(new.list=new.bf[[bf2]])
    new.bf[[bf2]] <- check.for.contiguous.OR.statments(new.list=new.bf[[bf2]])
    new.bf[b] <- NULL
  }
  return(new.bf)
}
