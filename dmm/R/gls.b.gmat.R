gls.b.gmat <-
function(am,vinv,dme.explist,vmatblock){
#  gls.b.gmat() - gls estimate of b using given Vinv matrix
#                 multiv2 version
#               - plus gmat which is like emat but based on M for gls$b
#V does not have to be positive definite ( uses ginv() )
#
   yblock <- matrix(am$y, am$n * am$l, 1) # yblock is l blocks each nx1
    xblock <- kronecker(diag(am$l), am$x) # xblock is lxl blocks each nxk
    vb <- ginv(t(xblock) %*% vinv %*% xblock)
    bblock <- vb %*% t(xblock) %*% vinv %*% yblock
#  unblock b to return
    b <- matrix(bblock, am$k,am$l)
#
    hblock <- xblock %*% vb %*% t(xblock) %*% vinv
    mblock <- diag(am$n * am$l) - hblock
#   yhatblock <- hblock %*% yblock
#   rblock <- mblock %*% yblock
#
    tmp <- matrix(0,am$n * am$l, am$n * am$l)
    gmat <- array(0,dim=c(am$n * am$n * am$l * am$l,am$v),
                    dimnames = list(NULL,dimnames(dme.explist$emat)[[2]]))
      for(iv in 1:am$v) {
        tmp <- mblock %*% matrix(vmatblock[,iv],am$n * am$l, am$n * am$l) %*% t(mblock)
#       gmat[ ,iv] <- as.vector(mblock %*% matrix(vmatblock[,iv],am$n * am$l, am$n * am$l) %*% t(mblock))  # prev version, gmat[,] in wrong order
        # reorder gmat to be block vectorized over 1st subscript
        gbeg <- 1
        for(il in 1 : am$l) {
          ibeg <- 1 + (il-1) * am$n
          iend <- il * am$n
          for(jl in 1: am$l) {
            jbeg <- 1 + (jl-1) * am$n
            jend <- jl * am$n
            gend <- gbeg + am$n * am$n - 1
            gmat[gbeg:gend,iv] <- tmp[ibeg:iend,jbeg:jend]
            gbeg <- gend + 1
          }
        }
      }
# 
    outlist <- list(b=b, bblock=bblock, vb=vb, yblock=yblock, xblock=xblock, gmat=gmat)
    return(outlist)
}
