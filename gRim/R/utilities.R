.rhsf2char <- function(f) {
  unlist(rhsf2list(f))
}

.list2pairs <- function(x){
  if (length(x)>0){
    if (!inherits(x,"list"))
      x <- list(x)
    x <- lapply(x, function(zzz) lapply(combn(zzz,2, simplify=FALSE), sort))
    x <- unlist(x, recursive=FALSE)
  }
  x
}

.toString <- function(x, col=' '){
  paste(x, collapse=col)
}

## Remove e and all higher order terms containing e from gen
##
## reduceSet(c(1,2,3,4),c(1,2))
## reduceSet(c(1,2,3,4),c(1,2,3))
.reduceSet <- function(gen,e){  
  if (length(e)==2){
    i <- match(e,gen)
    ans <- list(gen[-i[1]], gen[-i[2]])
  } else {
    le <- length(e)
    lx <- length(gen)
    x2 <- unlist(lapply(le:(lx-1), function(l) combn(gen,l, simplify=FALSE)),
                 recursive=FALSE)
    
    i <- isin(x2, e, TRUE)
    ans <- x2[i==0]
  }
  ans
}

## Remove term e from each element of glist
##
## delete.term(list(c(1,2,3,4),c(2,3,4,5)),c(1,2))
## delete.term(list(c(1,2,3,4),c(2,3,4,5)),c(1,2,3))

.delete.term <-  function(glist, e){
  idx<-isin(glist,e,TRUE)
  zzz <- unlist(
                lapply(1:length(glist), function(i){
                  if (idx[i]==1)
                    .reduceSet(glist[[i]],e)
                  else
                    glist[i]
                }), recursive=FALSE)
  ans <- removeRedundant(zzz)
  ans
}

## Add e interaction to x
##
.add.term <- function(glist,e){
  removeRedundant(c(glist,list(e)))
}
