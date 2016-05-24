summary.cna <- function(object, verbose=TRUE, ...) {

  ## summary.cna(net)
  ## y <- summary.cna(net, file="tmp.tbl", col.names=FALSE, append=T)
  ## or just write 'y' to file!
  
  if( !"cna" %in% class(object) ) {
    stop("Input should be a cna network object")
  }

  size <- table(object$communities$membership)
  id   <- names(size)
  memb <- sapply(id, function(i) { which(object$communities$membership==i) } )
  ## NOTE: Perhaps the memb should be names() of which inds
  ##       rather than the inds themselves as it is curently?
  ##memb <- sapply(id, function(i) { names(which(object$membership==i)) } )

  ##- Format as condensed vector for printing
  if( is.numeric(unlist(memb)) ) {
    members <- rep(NA, length(id))
    for(i in 1:length(id)) {
      b <- bounds(memb[[i]])[,c("start","end"),drop=FALSE]
      
      single.member <- NULL
      for(a in 1:dim(b)[1]){
        if(b[a,1] != b[a,2]){
          single.member[a] <- paste0(b[a,1], ":", b[a,2])
        }
        else{
          single.member[a] <- b[a,1]
        }
      }
      if(length(single.member)>1){
        members[i] <- paste0("c(", paste0(single.member, collapse=", "), ")")
      }
      else{
        members[i] <- single.member
      }
    }
  } else{
    ##- non numeric vectors can not be condensed
    members <- unlist(lapply(memb, paste, collapse=", "))
  }
  
  ## Output silently as a list  
  tbl <- data.frame( id=as.numeric(id), 
              size=as.numeric(size), 
              members=members,
              stringsAsFactors=FALSE )

  y <- list("id"=id, "size"=size, "members"=memb, "tbl"=tbl)
  
  if(verbose) { print.data.frame(tbl, row.names=FALSE) }

  return(y)
}

