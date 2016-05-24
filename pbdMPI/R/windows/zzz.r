### Lastest load into a package.

# .First.lib <- function(lib, pkg){
# } # End of .First.lib().

# .Last.lib <- function(libpath){
# } # End of .Last.lib().

.onLoad <- function(libname, pkgname){
  library.dynam("pbdMPI", pkgname, libname)

  # if(! is.loaded("spmd_initialize", PACKAGE = "pbdMPI")){
  #   stop("pbdMPI is not loaded.")
  # }

  # if(! exists(".__DISABLE_MPI_INIT__", envir = .GlobalEnv)){
  #   assign(".__DISABLE_MPI_INIT__", FALSE, envir = .GlobalEnv)
  # }

  if(exists(".__MPI_APTS__", envir = .GlobalEnv)){
    rm(list = c(".__MPI_APTS__"), envir = .GlobalEnv)
  }

  .Call("spmd_initialize", PACKAGE = "pbdMPI")

  ### For seed.
  # if(! exists(".lec.Random.seed.table", envir = .GlobalEnv) &&
  #    ! exists(".Random.seed", envir = .GlobalEnv)){
  if(! exists(".Random.seed", envir = .GlobalEnv)){
    # seed <- as.integer(Sys.getpid() + Sys.time())
    seed <- as.integer(runif(6, 1L, 2147483647L))
    seed <- .Call("spmd_bcast_integer", seed, 0L, 0L, PACKAGE = "pbdMPI")
    # seed <- rep(seed, 6)

    comm.size <- .Call("spmd_comm_size", 0L, PACKAGE = "pbdMPI")
    comm.rank <- .Call("spmd_comm_rank", 0L, PACKAGE = "pbdMPI")
    names <- as.character(0:(comm.size - 1))
    name <- as.character(comm.rank)

    invisible(eval(.lec.old.kind <- RNGkind(), envir = .GlobalEnv))
    invisible(eval(.lec.SetPackageSeed(seed), envir = .GlobalEnv))
    invisible(eval(.lec.CreateStream(names), envir = .GlobalEnv))
    invisible(eval(.lec.CurrentStream(name), envir = .GlobalEnv))
  }

  ### Preload to global environment.
  invisible(eval(parse(text = "pbdMPI:::.mpiopt_init()")))

  invisible()
} # End of .onLoad().

.onUnload <- function(libpath){
  pbdMPI::spmd.finalize(mpi.finalize = FALSE)
  library.dynam.unload("pbdMPI", libpath)

  invisible()
} # End of .onUnload().

.onAttach <- function(libname, pkgname){
#   if(interactive()){
#     pbdR.msg <-
# "
#     pbdR is mainly designed for batch environment or in SPMD programming.
#     Please utilize High Performance Computing Resource or consider an
#     allociation via ``http://www.nics.tennessee.edu/getting-an-allocation''
#     with XSEDE at ``http://www.xsede.org/overview''.
# "
#   } else{
#     pbdR.msg <-
# "
#     High Performance Computing Resource is available via
#     ``http://www.nics.tennessee.edu/getting-an-allocation''
#     with XSEDE at ``http://www.xsede.org/overview''.
# "
#   }
#   if(spmd.comm.rank(0L) == 0){
#     packageStartupMessage(pbdR.msg)
#   }
  invisible()
} # End of .onAttach().
