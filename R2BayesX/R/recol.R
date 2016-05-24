recol <-
function(x)
{
  xattr <- attributes(x)
  xattrn <- names(xattr)
  rn <- rownames(x)
  ok <- FALSE
  if(any(id <- grepl("10%", colnames(x)))) {
    ok <- TRUE
    cn <- colnames(x)
    cn <- cn[!id]
    x <- x[,!id]
    if(!is.matrix(x)) {
      x <- matrix(x, nrow = 1L)
      rownames(x) <- rn
      colnames(x) <- cn
    }
  }
  if(any(id <- grepl("90%", colnames(x)))) {
    ok <- TRUE
    cn <- colnames(x)
    cn <- cn[!id]
    x <- x[,!id]
    if(!is.matrix(x)) {
      x <- matrix(x, nrow = 1L)
      rownames(x) <- rn
      colnames(x) <- cn
    }
  }
  if(ok)
    for(i in 1L:length(xattrn))
      if(xattrn[i] != "dim" && xattrn[i] != "dimnames")
        attr(x, xattrn[i]) <- xattr[[i]]

  return(x)
}

