# Evaluates the splines of the indicated branch at the original projection index or.pi

lpc.spline.eval <- function(lpcsl,or.pi,branch=0){
                          
    if (length(branch)==1)
       branch <- rep(branch,length(branch))
    result <- matrix(nrow=length(or.pi), ncol=length(lpcsl[[1]]$splinefun))
    branches <- unique(branch)
    for (cur.branch in branches) {
      cur.or.pi <- or.pi[branch==cur.branch]
      for (j in 1:length(lpcsl[[cur.branch+1]]$splinefun))
         result[branch==cur.branch,j] <- lpcsl[[cur.branch+1]]$splinefun[[j]](cur.or.pi)  # here cur.op.pi instead of or.pi
      }
    result
  }
