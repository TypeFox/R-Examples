.subst <- function(x, i.val){
                                        #vv <- c("xyz[i+1]tyu", "xx[i]")
                                        #x <- c("xyz[i+1]tyu", "xx[i]","kkkk")
                                        #x <- c("xyztyu", "xx","kkkk")
  with.brack <- grep("\\[",x)
  vv <- x[with.brack]
  
  if (length(vv)>0){
    idx.vec <- gsub("[^\\[]*\\[([^\\]*)\\].*", "\\1", vv)
    idx.exp <- parse(text=idx.vec)
    idx.val <- unlist(lapply(idx.exp, eval, list(i=i.val)))
    vv2 <- list()
    for (ii in seq_along(idx.val)){
      vv2[[ii]] <- gsub("\\[([^\\]*)\\]", idx.val[ii], vv[ii])
    }
    vv2 <- unlist(vv2)
    x[with.brack] <- vv2
  }
  x
}

.do.one <- function(plist1, i.val){
  pp <- lapply(plist1, function(xx){
    xx$vpa <- .subst(xx$vpa,i.val)
    xx
  }) 
  pp 
}

repeatPattern <- function(plist, instances, unlist=TRUE){
  ans <- list()
  for (ii in seq_along(instances)){
    ans[[ii]] <- .do.one(plist, instances[[ii]])
  }
  if (unlist)
    ans <- unlist(ans, recursive=FALSE)

  ans
}
