layout_set <- function(mat, index)
  {
    mfg <- par("mfg")
    mfg[1:2] = which(mat==index, arr.ind=TRUE)[1,]
    par(mfg=mfg)
    invisible(mfg)
  }

layout_show <- function(mat)
  {
    graphics::layout.show( max(mat) )
    mat
  }
