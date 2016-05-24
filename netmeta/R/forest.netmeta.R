forest.netmeta <- function(x,
                           pooled=ifelse(x$comb.random, "random", "fixed"),
                           reference.group=x$reference.group,
                           leftcols="studlab",
                           leftlabs="Treatment",
                           smlab=NULL,
                           sortvar=x$seq,
                           ...){

  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")
  
  ipool <- charmatch(tolower(pooled), c("fixed", "random"), nomatch = NA)
  ##
  if (is.na(ipool)) 
        stop("Argument 'pooled' should be \"fixed\" or \"random\"")
  ##
  pooled <- c("fixed", "random")[ipool]
  
  if (pooled=="fixed"){
    TE   <- x$TE.fixed
    seTE <- x$seTE.fixed
    if (is.null(smlab))
      smlab <- "Fixed Effect Model"
  }
  ##
  if (pooled=="random"){
    TE   <- x$TE.random
    seTE <- x$seTE.random
    if (is.null(smlab))
      smlab <- "Random Effects Model"
  }

  labels <- colnames(TE)
  ##
  if (!is.null(sortvar)){
    if (is.character(sortvar)){
      seq <- setseq(sortvar, labels)
      TE <- TE[seq,seq]
      seTE <- seTE[seq,seq]
    }
    else{
      o <- order(sortvar)
      TE <- TE[o,o]
      seTE <- seTE[o,o]
    }
  }
  
  
  if (reference.group==""){
    warning("First treatment used as reference as argument 'reference.group' is unspecified.")
    reference.group <- labels[1]
  }
  else
    reference.group <- setref(reference.group, labels)
  ##
  TE.b <- TE[,colnames(TE)==reference.group]
  seTE.b <- seTE[,colnames(seTE)==reference.group]

  m1 <- metagen(TE.b, seTE.b, sm=x$sm,
                studlab=colnames(TE), warn=FALSE)
  
  forest(m1,
         comb.fixed=FALSE, comb.random=FALSE,
         hetstat=FALSE,
         leftcols=leftcols,
         leftlabs=leftlabs,
         smlab=smlab,
         ...)
  
  invisible(NULL)
}
