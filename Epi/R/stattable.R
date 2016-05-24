stat.table <- function(index, contents=count(), data, margins=FALSE)
{
  ## Get names of indices, using either user-supplied names, or by
  ## deparsing. Some of this code is stolen from the list.names()
  ## function in table()
  index.sub <- substitute(index)
  index <- if (missing(data))
    eval(index)
  else
    eval(index.sub, data)

  deparse.name <- function(x) if(is.symbol(x)) as.character(x) else ""
  if (is.list(index)) {
    if (is.call(index.sub)) {
      ## List constructed in the call to stat.table
      index.names <- names(index.sub)
      fixup <- if (is.null(index.names))
        seq(along = index.sub)
      else index.names == ""
      dep <- sapply(index.sub[fixup], deparse.name)
      if (is.null(index.names))
        index.labels <- dep
      else {
        index.labels <- index.names
        index.labels[fixup] <- dep
      }
      index.labels <- index.labels[-1]
    }
    else {
      ## List constructed outside the call to stat.table
      index.labels <- if (!is.null(names(index))) {
        names(index)
      }
      else {
        rep("", length(index))
      }
    }
  }
  else {
    ## A single vector
    index.labels <- deparse.name(index.sub)
  }
  
  ## Coerce index to list of factors
  if (!is.list(index))
    index <- list(index)
  index <- lapply(index, as.factor)

  ## Coerce contents to an expression representing a list of function calls
  contents <- substitute(contents)
  if (!identical(deparse(contents[[1]]), "list")) {  
    contents <- call("list", contents)
  }
  
  ## Check that functions in the contents list are valid
  valid.functions <- c("count","mean","weighted.mean","sum","quantile",
                       "median","IQR","max","min","ratio","percent","sd")
  table.fun <- character(length(contents) - 1)
  for (i in 2:length(contents)) {
    if (!is.call(contents[[i]]))
      stop("contents must be a list of function calls")
    FUN <- deparse(contents[[i]][[1]])
    if (!FUN %in% valid.functions)
      stop(paste("Function",FUN,"not permitted in stat.table"))
    else
      table.fun[i-1] <- FUN
  }
  
  ## Label the contents by default with the function call
  ## But if user has supplied a label use that instead.
  stat.labels <- sapply(contents, deparse)[-1]
  content.names <- names(contents)
  if (!is.null(content.names)) {
    for (i in 2:length(content.names)) {
      if (nchar(content.names[i]) > 0)
        stat.labels[i-1] <- content.names[i]
    }
  }

  ##Define the allowed tabulation functions
  
  count <- function(id){
      if (missing(id)) {
          id <- seq(along=index[[1]])
      }
      y <- tapply(id, INDEX=subindex, FUN=function(x) length(unique(x)))
      y[is.na(y)] <- 0
      return(y)
  }
  mean <- function(x, trim=0, na.rm=TRUE) {
    tapply(x, INDEX=subindex, FUN = base::mean, trim=trim, na.rm=na.rm)
  }
  weighted.mean <- function(x,w,na.rm=TRUE) {
    tapply(x, INDEX=subindex, FUN=stats::weighted.mean, w=w, na.rm=na.rm)
  }
  sum <- function(...,na.rm=TRUE) {
    tapply(..., INDEX=subindex, FUN = base::sum, na.rm=na.rm)
  }
  quantile <- function(x, probs, na.rm=TRUE,names=TRUE,type=7,...) {
    if (length(probs) > 1)
      stop("The quantile function only accepts scalar prob values within stat.table")
    tapply(x, INDEX=subindex, FUN = stats::quantile, probs=probs,
           na.rm=na.rm,names=names,type=type,...)
  }
  median <- function(x, na.rm=TRUE) {
    tapply(x, INDEX=subindex, FUN = stats::median, na.rm=na.rm)
  }
  IQR <- function(x, na.rm=TRUE) {
    tapply(x, INDEX=subindex, FUN= stats::IQR, na.rm=na.rm)
  }
  max <- function(..., na.rm=TRUE) {
    tapply(..., INDEX=subindex, FUN = base::max, na.rm=na.rm)
  }
  min <- function(..., na.rm=TRUE) {
    tapply(..., INDEX=subindex, FUN = base::min, na.rm=na.rm)
  }
  ratio <- function(d,y,scale=1, na.rm=TRUE) {
    if (length(scale) != 1)
      stop("Scale parameter must be a scalar")
    if (na.rm) {
      w <- (!is.na(d) & !is.na(y))
      tab1 <- tapply(d*w, INDEX=subindex, FUN=base::sum, na.rm=TRUE)
      tab2 <- tapply(y*w, INDEX=subindex, FUN=base::sum, na.rm=TRUE)
    }
    else {
      tab1 <- tapply(d, INDEX=subindex, FUN=base::sum, na.rm=FALSE)
      tab2 <- tapply(y, INDEX=subindex, FUN=base::sum, na.rm=FALSE)
    }
    return(scale*tab1/tab2)
  }
  percent <- function(...) {
    x <- list(...)
    if (length(x) == 0)
      stop("No variables to calculate percent")
    x <- lapply(x, as.factor)
    n <- count()
    ## Work out which indices to sweep out
    sweep.index <- logical(length(subindex))
    for (i in seq(along=subindex)) {
      sweep.index[i] <- !any(sapply(x,identical,subindex[[i]]))
    }
    if (!any(sweep.index)) {
      return(100*n/base::sum(n, na.rm=TRUE))
    }
    else {
      margin <- apply(n,which(sweep.index),base::sum, na.rm=TRUE)
      margin[margin==0] <- NA
      return(100*sweep(n, which(sweep.index), margin,"/"))
    }
  }
  sd <- function (..., na.rm = TRUE) 
  {
      tapply(..., INDEX=subindex, FUN = stats::sd, na.rm=na.rm)
  }

  ##Calculate dimension of the main table, excluding margins
  n.dim <- length(index)
  tab.dim <- sapply(index, nlevels)
               
  ##Sort out margins
  if (length(margins) == 1)
    margins <- rep(margins, n.dim)
  else if(length(margins) != n.dim)
    stop("Incorrect length for margins argument")

  ##Create grid of all possible subtables.
  fac.list <- vector("list", n.dim)
  for (i in 1:n.dim) {
    fac.list[[i]] <- if (margins[i]) c(0,1) else 1
  }
  subtable.grid <- as.matrix(expand.grid(fac.list))

  ##Fill in the subtables
  ans.dim <- c(length(contents)-1, tab.dim + margins)
  ans <- numeric(prod(ans.dim))
  for (i in 1:nrow(subtable.grid)) {
    ##in.subtable is a logical vector indicating which dimensions are
    ##in the subtable (i.e. which have not been marginalized out)
    in.subtable <- as.logical(subtable.grid[i,])
    llim <- rep(1,n.dim) + ifelse(in.subtable,rep(0,n.dim),tab.dim)
    ulim <- tab.dim + ifelse(in.subtable,rep(0,n.dim),rep(1, n.dim))
    subindex <- index[in.subtable]
    if (length(subindex) == 0) {
        ## Marginalizing out all dimensions
        subindex <- list(rep(1, length(index[[1]])))
    }
    subtable.list <- if(missing(data))
        ##eval(contents, parent.frame())
        eval(contents)
    else
      eval(as.expression(contents), data)
    for (j in 1:length(subtable.list)) {
      ans[array.subset(ans.dim,c(j,llim),c(j,ulim))] <- subtable.list[[j]]
    }
  }
  
  ans <- array(ans, dim=ans.dim)
  ans.dimnames <- lapply(index, levels)
  names(ans.dimnames) <- index.labels
  for (i in 1:length(index)) {
    if (margins[i])
      ans.dimnames[[i]] <- c(ans.dimnames[[i]], "Total")
  }
  dimnames(ans) <- c(list("contents"=stat.labels), ans.dimnames)
  attr(ans, "table.fun") <- table.fun
  class(ans) <- c("stat.table", class(ans))
  return(ans)
}


array.subset <- function(dim,lower,upper)
{
  ##Returns a logical array of dimension dim for which elements in the range
  ##[lower[1]:upper[1], lower[2]:upper[2],...] are TRUE and others FALSE

  ##Check validity of arguments (but assume everything is an integer)
  ndim <- length(dim)
  if (length(lower) != ndim || length(upper) != ndim) {
    stop("Length of lower and upper limits must match dimension")
  }
  if (any(lower > upper) || any(lower < 1) || any(upper > dim)) {
    stop("Invalid limits")
  }

  ##The math is easier if we index arrays from 0 rather than 1
  lower <- lower - 1
  upper <- upper - 1

  N <- prod(dim)
  ans <- rep(TRUE, N)
  for (i in 1:N) {
    l <- i - 1
    for (d in 1:ndim) {
      k <- l %% dim[d] #k is the index of the ith element in dimension d
      if (k < lower[d] || k > upper[d]) {
        ans[i] <- FALSE
        break
      }
      l <- l %/% dim[d] 
    }
  }
  return(array(ans, dim))
}

split.to.width <- function(x,width)
{
  ## Splits a string into a vector so that each element has at most width
  ## characters.  If width is smaller than the length of the shortest word
  ## then the latter is used instead
  
  x.split <- strsplit(x,split=" ")[[1]]
  width <- max(c(width,nchar(x.split)))
  y <- character(0)
  imin <- 1
  n <- length(x.split)
  for (i in 1:n) {
    cum.width <- if(i==n) {
      Inf
    }
    else {
      sum(nchar(x.split[imin:(i+1)])) + (i - imin + 1)
    }
    if (cum.width > width) {
      y <- c(y,paste(x.split[imin:i], collapse=" "))
      imin <- i + 1
    }
  }
  return(y)
}

pretty.print.stattable.1d <- function(x, width, digits)
{
  ##Pretty printing of 1-D stat.table

  if (length(dim(x)) != 2)
    stop("Cannot print stat.table")
  ncol <- nrow(x)
  
  col.width <- numeric(ncol+1)
  col.header <- vector("list",ncol+1)
  n.header <- integer(ncol+1)
  print.list <- vector("list",ncol+1)
  
  ##First column
  col.header[[1]] <- split.to.width(names(dimnames(x))[2], width)
  n.header[1] <- length(col.header[[1]])
  col1 <- format(c(col.header[[1]],dimnames(x)[[2]]), justify="left")
  col.header[[1]] <- col1[1:n.header[1]]
  print.list[[1]] <- col1[-(1:n.header[1])]
  col.width[1] <- nchar(col.header[[1]][1])
  
  ##Other columns
  for (i in 2:(ncol+1)) {
    col.header[[i]] <- split.to.width(dimnames(x)[[1]][i-1], width)
    n.header[i] <- length(col.header[[i]])
    this.col  <- formatC(x[i-1,],width=width,
                         digits=digits[attr(x,"table.fun")[i-1]], "f")
    this.col <- format(c(col.header[[i]],this.col),justify="right")
    col.width[i] <- nchar(this.col[1])
    col.header[[i]] <- this.col[1:n.header[i]]
    print.list[[i]] <- this.col[-(1:n.header[i])]
  }
    
  ##
  table.width <- sum(col.width) + ncol + 3
  max.n.header <- max(n.header)
  cat(" ",rep("-",table.width)," \n",sep="")
  for(i in 1:max.n.header) {
    cat(" ")
    for(j in 1:length(print.list)) {
      if (i <= n.header[j]) {
        cat(col.header[[j]][i])
      }
      else {
        cat(rep(" ", col.width[[j]]),sep="")
      }
      if (j==1)
        cat("   ")
      else
        cat(" ")
    }
    cat(" \n")
  }
  cat(" ",rep("-",table.width)," \n",sep="")

  for (i in 1:length(print.list[[1]])) {
    cat(" ")
    if (pmatch("Total",print.list[[1]][i],nomatch=0)) {
      ##Add a blank line before the total
      cat(rep(" ",col.width[1]+1),"  ",rep(" ",sum(col.width[-1])+ncol),
          " \n ",sep="")
    }
    for (j in 1:length(print.list)) {
      cat(print.list[[j]][i])
      if (j == 1) {
        cat("   " )
      }
      else {
        cat(" ")
      }
    }
    cat(" \n")
  }
  cat(" ",rep("-",table.width)," \n",sep="")
  return(invisible(x))
}

pretty.print.stattable.2d <- function(x, width, digits)
{
  ##Pretty printing of 2-Dimensional stat.table

  if (length(dim(x)) != 3)
    stop("Cannot print stat.table")
  nstat <- dim(x)[1]
  ncol <- dim(x)[3]
  nrow <- dim(x)[2]
  
  col.width <- numeric(ncol+1)
  col.header <- vector("list",ncol+1)
  n.header <- integer(ncol+1)
  print.list <- vector("list",ncol+1)

  ##First column
  col.header[[1]] <- split.to.width(names(dimnames(x))[2], width)
  n.header[1] <- length(col.header[[1]])
  col1 <- format(c(col.header[[1]],dimnames(x)[[2]]), justify="left")
  col.header[[1]] <- col1[1:n.header[1]]
  print.list[[1]] <- col1[-(1:n.header[1])]
  col.width[1] <- nchar(col.header[[1]][1])

  ##Other columns
  for (i in 2:(ncol+1)) {
    col.header[[i]] <- split.to.width(dimnames(x)[[3]][i-1], width)
    n.header[i] <- length(col.header[[i]])
    this.col <- matrix("", nrow=nstat,ncol=nrow)
    for (j in 1:nstat) {
      z <- x[j,,i-1]
      this.col[j,]  <- formatC(z, width=width, format="f", 
                               digits=digits[attr(x,"table.fun")[j]])
      ##      this.col[j,] <- formatC(z, width=width, digits=digits,
      ##                       format=ifelse(identical(round(z),z),"d","f"))
    }
    this.col <- format(c(col.header[[i]],this.col),justify="right")
    col.width[i] <- nchar(this.col[1])
    col.header[[i]] <- this.col[1:n.header[i]]
    print.list[[i]] <- this.col[-(1:n.header[i])]
  }

  ##Correct first column for multiple stats
  if (nstat > 1) {
    pl1 <- print.list[[1]]
    print.list[[1]] <- rep(paste(rep(" ",col.width[1]),collapse=""),nstat*nrow)
    print.list[[1]][1 + nstat*((1:nrow)-1)] <- pl1
  }

  table.width <- sum(col.width) + ncol + 3
  max.n.header <- max(n.header)
  cat(" ",rep("-",table.width)," \n",sep="")
  ## Supercolumn header
  super.header <- names(dimnames(x))[3]
  npad <- sum(col.width[-1]) + ncol + 1 - nchar(super.header)
  if (npad >= 0) {
    cat(" ",rep(" ",col.width[1]),"  ",sep="")
    cat(rep("-",floor(npad/2)),sep="")
    cat(super.header)
    cat(rep("-",ceiling(npad/2))," \n",sep="")
  }
  ## Headers
  for(i in 1:max.n.header) {
    cat(" ")
    for(j in 1:length(print.list)) {
      if (i <= n.header[j]) {
        cat(col.header[[j]][i])
      }
      else {
        cat(rep(" ", col.width[[j]]),sep="")
      }
      if (j==1)
        cat("   ")
      else
        cat(" ")
    }
    cat(" \n")
  }
  cat(" ",rep("-",table.width)," \n",sep="")
  ## Body of table
  blank.line <- function() {
    cat(" ",rep(" ",col.width[1]+1),"  ",rep(" ",sum(col.width[-1])+ncol),
        " \n",sep="")
  }
  for (i in 1:length(print.list[[1]])) {
    if (pmatch("Total",print.list[[1]][i],nomatch=0)) {
      ##Add a blank line before the total
      blank.line()
    }
    cat(" ")
    for (j in 1:length(print.list)) {
      cat(print.list[[j]][i])
      if (j == 1) {
        cat("   " )
      }
      else {
        cat(" ")
      }
    }
    cat(" \n")
    if (nstat > 1 && i %% nstat == 0 && i != length(print.list[[1]])) {
      ##Separate interleaved stats
      blank.line()
    }
  }
  cat(" ",rep("-",table.width)," \n",sep="")
  return(invisible(x))
}

print.stat.table <- function(x, width=7,digits,...)
{
  fun.digits <- c("count"=0,"mean"=2,"weighted.mean"=2,"sum"=2,"quantile"=2,
                  "median"=2,"IQR"=2,"max"=2,"min"=2,"ratio"=2,"percent"=1,
                  "sd"=2)
  if (!missing(digits)) {
    if (is.null(names(digits))) {
      if (length(digits) > 1)
        stop("digits must be a scalar or named vector")
      else
        fun.digits[1:length(fun.digits)] <- digits
    }
    else {
      fun.digits[names(digits)] <- digits
    }
  }
  if (length(dim(x)) == 2)
    pretty.print.stattable.1d(x, width, fun.digits)
  else if (length(dim(x)) == 3)
    pretty.print.stattable.2d(x, width, fun.digits)
  else
    NextMethod("print",...)
}

## Satisfy QA checks by defining  these functions. But if we never
## export them they can't be used directly.

count <- function(id)
{
}

ratio <- function(d, y, scale=1, na.rm=TRUE)
{
}

percent <- function(...)
{
}
