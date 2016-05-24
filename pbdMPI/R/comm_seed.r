### Seed functions for random number generators.

comm.set.seed <- function(seed, diff = FALSE, state = NULL,
    comm = .pbd_env$SPMD.CT$comm){
  if(exists(".lec.Random.seed.table", envir = .GlobalEnv)){
    comm.end.seed(comm)
  }
  if(length(seed) == 1){
    seed <- rep(seed, 6)
  }
  seed <- as.integer(seed[1:6])
  seed <- spmd.bcast.integer(seed, rank.source = 0L, comm = comm)

  if(diff){
    names <- as.character(0:(comm.size(comm) - 1))
    name <- as.character(comm.rank(comm))
  } else{
    names <- "0"
    name <- "0"
  }

  invisible(eval(.lec.old.kind <- RNGkind(), envir = .GlobalEnv))
  invisible(eval(.lec.SetPackageSeed(seed), envir = .GlobalEnv))
  invisible(eval(.lec.CreateStream(names), envir = .GlobalEnv))
  if(! is.null(state)){
    invisible(eval(.lec.SetSeed(name, state), envir = .GlobalEnv))
  }
  invisible(eval(.lec.CurrentStream(name), envir = .GlobalEnv))
  invisible()
} # End of comm.set.seed().

comm.end.seed <- function(comm = .pbd_env$SPMD.CT$comm){
  name <- get(".lec.Random.seed.table", envir = .GlobalEnv)$name
  old.kind <- try(get(".lec.old.kind", envir = .GlobalEnv), silent = TRUE)
  if(class(old.kind) == "try-error"){
    old.kind <- RNGkind()
  }

  if(! is.null(name)){
    invisible(eval(.lec.CurrentStreamEnd(old.kind), envir = .GlobalEnv))
    invisible(eval(.lec.DeleteStream(name), envir = .GlobalEnv))
    # rm.list <- c(".lec.Random.seed.table", ".Random.seed")
    # invisible(eval(rm(list = rm.list, envir = .GlobalEnv)))
  }

  invisible()
} # End of comm.end.seed().

comm.reset.seed <- function(comm = .pbd_env$SPMD.CT$comm){
  seed.table <- get(".lec.Random.seed.table", envir = .GlobalEnv)
  if(is.null(seed.table)){
    comm.stop("seed.table is not found.", comm = comm)
  }

  old.kind <- try(get(".lec.old.kind", envir = .GlobalEnv), silent = TRUE)
  if(class(old.kind) == "try-error"){
    old.kind <- RNGkind()
  }

  if(length(seed.table$name) == comm.size(comm)){
    name <- seed.table$name[comm.rank(comm) + 1]
  } else if(length(seed.table$name) == 1){
    name <- seed.table$name[1]
  } else{
    comm.stop("seed.table is incorrect.", comm = comm)
  }

  invisible(eval(.lec.CurrentStreamEnd(old.kind), envir = .GlobalEnv))
  invisible(eval(.lec.ResetStartStream(name), envir = .GlobalEnv))
  invisible(eval(.lec.CurrentStream(name), envir = .GlobalEnv))
  invisible()
} # End of comm.reset.seed().

comm.seed.state <- function(comm = .pbd_env$SPMD.CT$comm){
  seed.table <- get(".lec.Random.seed.table", envir = .GlobalEnv)
  if(is.null(seed.table)){
    comm.stop("seed.table is not found.", comm = comm)
  }

  old.kind <- try(get(".lec.old.kind", envir = .GlobalEnv), silent = TRUE)
  if(class(old.kind) == "try-error"){
    old.kind <- RNGkind()
  }

  if(length(seed.table$name) == comm.size(comm)){
    name <- seed.table$name[comm.rank(comm) + 1]
  } else if(length(seed.table$name) == 1){
    name <- seed.table$name[1]
  } else{
    comm.stop("seed.table is incorrect.", comm = comm)
  }

  ### Terminate the stream first to get the current state.
  invisible(.lec.CurrentStreamEnd(old.kind))
  ret <- .lec.GetState(name)

  ### Set the state back to stream (pretend as nothing happens).
  invisible(.lec.SetSeed(name, ret))
  invisible(.lec.CurrentStream(name))
  ret
} # End of comm.seed.state().

