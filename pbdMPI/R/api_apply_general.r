### These functions are supposed to run in SPMD, even when pbd.model = "mw".

build.index.matrix <- function(dim.X, MARGIN){
  tmp.dim <- dim.X[MARGIN]
  ret <- cbind(1:tmp.dim[1])
  if(length(tmp.dim) > 1){
    for(i in 2:length(tmp.dim)){
      if(tmp.dim[i] > 1){
        tmp <- ret
        for(j in 2:tmp.dim[i]){
          ret <- rbind(ret, tmp)
        }
      }
      ret <- cbind(ret, rep(1:tmp.dim[i], each = nrow(ret) / tmp.dim[i]))
    }
  }

  ret
} # End of build.index.matrix().

build.index.list <- function(dim.X, MARGIN){
  ret <- build.index.matrix(dim.X, MARGIN)
  ret <- lapply(1:nrow(ret), function(i) ret[i,])
  ret
} # End of build.index.matrix().

subset.by.index <- function(i.id, X, dim.X, MARGIN){
  tl <- length(dim.X)
  text.X <- "X["
  for(i in 1:tl){
    if(i %in% MARGIN){
      text.X <- paste(text.X, "i.id[MARGIN == ",  i, "]", sep = "")
    }
    if(i < tl){
      text.X <- paste(text.X, ",", sep = "")
    }
  }
  text.X <- paste(text.X, "]", sep = "")
  ret <- eval(parse(text = text.X))

  ret
} # End of subset.by.index().

array.to.list <- function(jid, X, dim.X, MARGIN){
  tmp.id <- build.index.list(dim.X, MARGIN)
  ret <- lapply(tmp.id[jid], subset.by.index, X, dim.X, MARGIN) 

  ret
} # End of array.to.list().

pbdApply.general <- function(X, MARGIN, FUN, ...,
    pbd.mode = c("mw", "spmd", "dist"),
    rank.source = .pbd_env$SPMD.CT$rank.root, comm = .pbd_env$SPMD.CT$comm){
  COMM.SIZE <- spmd.comm.size(comm)
  COMM.RANK <- spmd.comm.rank(comm)

  ### Redistributing data for MW.
  if(pbd.mode[1] == "mw"){
    MARGIN <- as.integer(MARGIN)
    new.X <- NULL 
    if(COMM.RANK == rank.source){
      if(! is.array(X)){
        stop("X should be an array")
      }
      alljid <- get.jid(prod(dim(X)[MARGIN]), comm = comm, all = TRUE)
      new.X <- lapply(alljid, array.to.list, X, dim(X), MARGIN) 
    }

    if(length(new.X) < COMM.SIZE){
      new.X <- c(new.X, rep(list(NULL), COMM.SIZE - length(new.X)))
    }

    new.X <- spmd.scatter.object(new.X, rank.source = rank.source, comm = comm)

    ### For dots.
    arg.dots <- list(...)
    if(COMM.RANK == rank.source){
      arg.dots <- rep(list(arg.dots), COMM.SIZE)
    }
    arg.dots <- spmd.scatter.object(arg.dots, rank.source = rank.source,
                                    comm = comm)
    if(length(arg.dots) > 0){
      names <- names(arg.dots)
      for(i in 1:length(arg.dots)){
        assign(names[i], arg.dots[[i]])
      }
    }
  } else if(pbd.mode[1] == "spmd"){
    alljid <- get.jid(prod(dim(X)[MARGIN]), comm = comm)
    new.X <- lapply(list(alljid), array.to.list, X, dim(X), MARGIN) 
    new.X <- unlist(new.X, recursive = FALSE)
  } else if(pbd.mode[1] == "dist"){
    new.X <- X
  } else{
    comm.stop("pbd.mode is not found.")
  }

  ### Run as SPMD.
  if(pbd.mode[1] != "dist"){
    ret <- lapply(new.X, FUN, ...)
  } else{
    ret <- apply(new.X, MARGIN, FUN, ...)
  }

  ### Gather data for MW.
  if(pbd.mode[1] == "mw"){
    ret <- spmd.gather.object(ret, rank.dest = rank.source, comm = comm)
    if(COMM.RANK != rank.source){
      ret <- NULL
    } else{
      ret <- unlist(ret, use.names = FALSE)
      dim(ret) <- dim(X)[MARGIN]
      dimnames(ret) <- dimnames(X)[MARGIN]
    }
  }

  ret
} # End of pbdApply.general().

