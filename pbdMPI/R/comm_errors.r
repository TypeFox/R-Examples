### Similar to stop(), warning(), and warnings().

comm.stop <- function(..., call. = TRUE, domain = NULL,
    all.rank = .pbd_env$SPMD.CT$print.all.rank,
    rank.print = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm,
    mpi.finalize = .pbd_env$SPMD.CT$mpi.finalize,
    quit = .pbd_env$SPMD.CT$quit){
  COMM.RANK <- spmd.comm.rank(comm)
  spmd.finalize(mpi.finalize = mpi.finalize)

  if(COMM.RANK %in% rank.print || all.rank == TRUE){
    args <- list(...)
    if(length(args) == 1L && inherits(args[[1L]], "condition")){
      cond <- args[[1L]]
      if(nargs() > 1L){
        warning("additional arguments ignored in comm.stop()")
      }
      message <- conditionMessage(cond)
      call <- conditionCall(cond)
      ### This is not allowed by "R CMD check".
      # .Internal(.signalCondition(cond, message, call))
      # .Internal(.dfltStop(message, call))

      ### WCC: I don't know what this if() for. The next should work,
      ###      but may not fit the original purpose.
      .External("api_R_stop", call, as.logical(call.), message,
                PACKAGE = "pbdMPI")
    } else{
      ### This is not allowed by "R CMD check".
      # .Internal(stop(as.logical(call.), .makeMessage(..., 
      #           domain = domain)))

      .External("api_R_stop", sys.call(-1), as.logical(call.),
                .makeMessage(..., domain = domain),
                PACKAGE = "pbdMPI")
    }
  }

  if(quit){
    q("no", status = 1, runLast = TRUE)
  }

  invisible()
} # End of comm.stop().

comm.warning <- function(..., call. = TRUE, immediate. = FALSE, domain = NULL,
    all.rank = .pbd_env$SPMD.CT$print.all.rank,
    rank.print = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  if(spmd.comm.rank(comm) %in% rank.print || all.rank == TRUE){
    args <- list(...)
    if(length(args) == 1L && inherits(args[[1L]], "condition")){
      cond <- args[[1L]]
      if(nargs() > 1L) 
        cat(gettext("additional arguments ignored in comm.warning()"), 
            "\n", sep = "", file = stderr())
        message <- conditionMessage(cond)
        call <- conditionCall(cond)
        withRestarts({
            ### This is not allowed by "R CMD check".
            # .Internal(.signalCondition(cond, message, call))
            # .Internal(.dfltWarn(message, call))

            ### WCC: I don't know what this if() for. The next should work,
            ###      but may not fit the original purpose.
            .External("api_R_warning", call, as.logical(call.),
                      as.logical(immediate.), message,
                      PACKAGE = "pbdMPI")
        }, muffleWarning = function() NULL)
        invisible(message)
    } else{
      ### This is not allowed by "R CMD check".
      # .Internal(warning(as.logical(call.), as.logical(immediate.), 
      #           .makeMessage(..., domain = domain)))

      ret <- .External("api_R_warning", sys.call(-1), as.logical(call.),
                       as.logical(immediate.),
                       .makeMessage(..., domain = domain),
                       PACKAGE = "pbdMPI")
      invisible(ret)
    }
  }
} # End of comm.warning().

comm.warnings <- function(...,
    all.rank = .pbd_env$SPMD.CT$print.all.rank,
    rank.print = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  if(spmd.comm.rank(comm) %in% rank.print || all.rank == TRUE){
    if (!exists("last.warning", envir = baseenv())) 
        return()
    last.warning <- get("last.warning", envir = baseenv())
    if (!(n <- length(last.warning))) 
        return()
    structure(last.warning, dots = list(...), class = "warnings")
  }
} # End of comm.warnings().

