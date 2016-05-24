## This file (R/truthTable.R) is part of QCA3 package
## copyright: HUANG Ronggui 2008-2011

findNoCSA <- function(x,y, noCSA=TRUE){
    nx <- length(x$solutions)
    ny <- length(y$solutions)
    idx <- expand.grid(x=seq(1,nx),y=seq(1,ny))
    ncsa <- apply(idx,1,FUN=function(ii) {
        csa <- CSA(x[ii[1]],y[ii[2]])$solutions
        if ((nc <-length(csa)) == 1) {
        if (is.null(csa[[1]])) nc <- 0
       }
     nc
    }
                  )
    ans <- cbind(idx,ncsa=ncsa)
    if (noCSA) ans <- ans[ans$ncsa==0,]
    ans
}

##bestSubsets(ans1,ans0)


findParsimonious <- function(x) {
    ans <- sapply(x$solutions, function(i) sum(!is.dontcare(i)))
    idx <- which(ans == min(ans))
    idx
}
