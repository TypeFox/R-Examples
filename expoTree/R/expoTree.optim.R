expoTree.optim <- function(forest,
      fix=rep(5,F),fix.val=rep(5,0),
      pars=vector(length=sum(!fix))*1,
      lo=rep(sum(!fix),0),hi=rep(sum(!fix),0),
      survival=TRUE,vflag=0,method="Nelder-Mead",control=list(trace=0)) 
{
  if (! is.list(forest)) stop("Must supply a list of trees.")

  fn <- function(par,forest,fix,fix.val,survival,vflag) {
    x <- fix.val
    cnt <- 1
    for (i in 1:5) {
      if (! fix[i]) {
        x[i] <- par[cnt]
        cnt <- cnt + 1
      }
    }
    if (x[3] < 0 & fix[3]) x[3] <- x[4]*(1.0/(-x[3])-1.0)
    print(par)
    print(cbind(lo,hi))
    if (any(par < lo | par > hi)) return(-Inf)
    lik <- sapply(forest,function(tree) {
      l <- -Inf
      if (x[1] > 0) {
        l <- runExpoTree(x,times=tree[,1],ttypes=tree[,2],survival=survival,vflag=vflag)
      } else {
        l <- infExpoTree(x[2:5],times=tree[,1],ttypes=tree[,2],survival=survival,vflag=vflag)
      }
      return(l)
    })
    return(sum(lik))
  }

  control <- list(fnscale=-1,parscale=pars)
  optim(par=pars,fn=fn,control=control,method=method,
        forest=forest,fix=fix,fix.val=fix.val,survival=survival,vflag=vflag)
}

