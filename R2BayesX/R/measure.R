measure <-
function(a, b, model = NULL, fixed = FALSE)
{
  err <- list()
  ## compare fixed effects
  if(!is.null(a[[1L]]$fixed.effects) && !is.null(b[[1L]]$fixed.effects)) {
    fa <- a[[1L]]$fixed.effects
    fb <- a[[1L]]$fixed.effects
    spez <- paste(c(model, "fixed.effects"), sep = "", collapse = "_")
    eval(parse(text=paste("err$", spez, "<- sum((fa-fb)^2)", sep = "")))
  } else return(NULL)
  if(!fixed) {
    ## compare smooth effects
    a <- a[[1L]]$effects
    b <- b[[1L]]$effects
    n <- length(a)
    m <- length(b)
		k <- 1L
    for(i in 1L:n) {
      A <- a[[i]]
      c1 <- colnames(A)
      la <- attr(A, "specs")$label
      A <- A[!duplicated(A[,1L]),]
      A <- A[order(A[,1L]),]
      for(j in 1L:m) {
        B <- b[[j]]
        c2 <- colnames(B)
        lb <- attr(B, "specs")$label
        B <- B[!duplicated(B[,1L]),]
        B <- B[order(B[,1L]),]
        if(all(c1[1L:2L] == c2[1L:2L]) && la==lb) {
          if((nrow(A) == nrow(B)) && (A[,1L]==B[,1L])) {
            start <- 2L
            if(length(c2 ) == 10L)
							start <- 3L
						if(length(c2) == 11L)
							start <- 4L
						for(d in start:(ncol(A) - 1L)) {
              spez <- c(model,c1[1L:(start - 1L)])
              spez <- paste(c(spez,c1[d]), sep = "", collapse = "_")
              eval(parse(text = paste("err$", spez, "<- sum((A[,d]-B[,d])^2)/nrow(A)", sep = "")))
            }
          } else stop("missmatching matrices!")
        }
      }
    }
  }
  if(sum(unlist(err)) > 0)
    cat("Detecting differences in ",model,": ",sum(unlist(err)),"\n",sep="")

  return(err)
}

