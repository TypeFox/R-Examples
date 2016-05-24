suppressPackageStartupMessages(library(pbdMPI, quietly=TRUE))
init()

if (comm.size() < 2) comm.stop("You need at least 2 MPI ranks for this test")

if (comm.rank() == 0){
  x <- TRUE
} else {
  x <- FALSE
}

y <- bcast(x)

comm.print(y, rank.print=1)

finalize()
