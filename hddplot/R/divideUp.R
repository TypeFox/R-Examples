"divideUp" <-
function(cl=rep(1:3, c(7,4,8)), nset=2, seed=NULL, balanced=TRUE){
    if(!is.null(seed))set.seed(seed)
    if(balanced){
      ord <- order(cl)
      ordcl <- cl[ord]
      gp0 <- rep(sample(1:nset), length.out=length(cl))
      gp <- unlist(split(gp0,ordcl), function(x)sample(x))
      gp[ord] <- gp
    } else
    gp <- sample(rep(1:nset, length.out=length(cl)))
    as.vector(gp)
  }

