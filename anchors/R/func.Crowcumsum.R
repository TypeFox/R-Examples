Crowcumsum <- function(x,n.row,n.col,n.block=1) {
  if(all(is.finite(x))) {
    return(.C("vecrowcumsum",x=as.double(x),as.integer(n.block),
              as.integer(n.col),as.integer(n.row),
              PACKAGE = "anchors")$x)
  } else {
    return(NULL)
  }
}
