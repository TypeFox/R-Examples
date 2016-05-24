### These functions are supposed to run in SPMD, even when pbd.model = "mw".

array.to.list.RC <- function(jid, X, dim.X, MARGIN){
  tl <- length(dim.X)
  text.X <- "X["
  for(i in 1:tl){
    if(i == MARGIN){
      text.X <- paste(text.X, "c(jid)", sep = "")
    }
    if(i < tl){
      text.X <- paste(text.X, ",", sep = "")
    }
  }
  text.X <- paste(text.X, "]", sep = "")
  ret <- eval(parse(text = text.X))

  dim.X[MARGIN] <- length(jid)
  dim(ret) <- dim.X
  ret
} # End of array.to.list.RC().

pbdApply.RC <- function(X, MARGIN, FUN, ...,
    pbd.mode = c("mw", "spmd", "dist"),
    rank.source = .pbd_env$SPMD.CT$rank.root, comm = .pbd_env$SPMD.CT$comm){
  MARGIN <- MARGIN[1]

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
      alljid <- get.jid(dim(X)[MARGIN], comm = comm, all = TRUE)
      new.X <- lapply(alljid, array.to.list.RC, X, dim(X), MARGIN) 
    }

    if(length(new.X) < COMM.SIZE){
      new.X <- c(new.X, rep(list(NULL), COMM.SIZE - length(new.X)))
    }

    new.X <- spmd.scatter.array(new.X, rank.source = rank.source, comm = comm)

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
    alljid <- get.jid(dim(X)[MARGIN], comm = comm)
    new.X <- sapply(alljid, array.to.list.RC, X, dim(X), MARGIN) 
  } else if(pbd.mode[1] == "dist"){
    new.X <- X
  } else{
    comm.stop("pbd.mode is not found.")
  }

  ### Run as SPMD.
  ret <- apply(new.X, MARGIN, FUN, ...)

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
} # End of pbdApply.RC().

