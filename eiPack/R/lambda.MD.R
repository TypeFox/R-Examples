lambda.MD <- function(object, columns, ret.mcmc = TRUE){

  get2 <- function(x) x[2]
  if (class(object) != "eiMD")
    stop("'object' must be output from 'ei.MD.bayes'")
  if (missing(columns) | length(columns) < 2)
    stop("'columns' requires at least two column names")
  cc <- object$draws$Cell.counts
  if (is.mcmc(cc)) { 
    tnames <- strsplit(colnames(cc), "ccount.")
    idx <- strsplit(sapply(tnames, get2), ".", fixed = TRUE)
    idx <- as.list(as.data.frame(matrix(unlist(idx), byrow = TRUE,
                                        nrow = length(idx), ncol = 
                                        length(idx[[1]]))))
    idx <- lapply(idx, as.character)
    idx <- lapply(idx, unique)
    sims <- nrow(cc)
    mcpar.cc <- mcpar(cc)
    cc <- array(t(cc), c(sapply(idx, length), sims),
                dimnames = list(idx[[1]], idx[[2]], 1:sims))
  }
  else { 
    idx <- names(object$draws$Cell.counts)[1:2]
    sims <- dim(cc)[3]
  }
  names(idx) <- c("rows", "columns")
  NG <- length(idx[[1]])
  NP <- length(idx[[2]])
  NI <- length(columns)

  lambda.out <- array(NA, c(NG, NI, sims),
                      dimnames = list(idx[[1]], columns, 1:sims))
  for (ii in idx[[1]]) {
    lambda.out[ii,,] <- t(t(cc[ii, columns,]) /
                          apply(cc[ii, columns,], 2, sum))
  }
  if (ret.mcmc){
    lambda.out <- t(matrix(lambda.out , NG*NI, sims))
    colnames(lambda.out) <-  paste("lambda", matrix(rep(idx[[1]],
                                                        NI),NG,NI),
                                   matrix(rep(columns, NG),NG,NI,
                                          byrow=TRUE) ,sep=".")
    lambda.out <- coda::mcmc(lambda.out)
    attr(lambda.out, "mcpar") <- mcpar.cc
  }

  class(lambda.out) <- c("lambdaMD", class(lambda.out))
  lambda.out
}


