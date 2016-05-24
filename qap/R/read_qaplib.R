
read_qaplib <- function(file) {
  if(!file.exists(file)) stop("file ", file, " does not exist!")

  dat <- as.integer(scan(file, quiet = TRUE))
  n <- dat[1]

  A <- matrix(dat[2:(n*n+1L)], ncol = n, nrow = n, byrow = TRUE)
  B <- matrix(dat[(n*n+2L):(n*n+2L+n*n-1L)], ncol = n, nrow = n, byrow = TRUE)

  # read solution if available
  sol <- NULL
  opt <- NULL
  file_sol <- paste(sub("^([^.]*).*", "\\1", file), ".sln", sep ="")
  if(file.exists(file_sol)) {
    dat <- scan(file_sol, quiet = TRUE)
    sol <- dat[-(1:2)]
    opt <- dat[2]
  }


  list(A=A, B=B, solution = sol, opt = opt)
}
