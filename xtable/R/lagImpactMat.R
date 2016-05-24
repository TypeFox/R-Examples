### This function is a copy of spdep:::lagImpactMat
### It has been copied because lagImpactMat is not exported by spdep
### There is no help available for lagImpactMat in spdep,
### so I have not provided any help, and I am unable to trace the author
###
lagImpactMat <- function (x, reportQ = NULL) 
{
  if (is.null(x$res)) {
    direct <- x$direct
    indirect <- x$indirect
    total <- x$total
  } else {
    direct <- x$res$direct
    indirect <- x$res$indirect
    total <- x$res$total
  }
  mat <- cbind(direct, indirect, total)
  colnames(mat) <- c("Direct", "Indirect", "Total")
  rownames(mat) <- attr(x, "bnames")
  if (!is.null(reportQ) && reportQ) {
    if (is.null(x$res)) {
      Qobj <- attr(x, "Qres")
    } else {
      Qobj <- attr(x$res, "Qres")
    }
    if (is.null(Qobj)) {
      warning("No impact components to report")
    } else {
      if (length(attr(x, "bnames")) == 1L) {
        Qobj$direct <- matrix(Qobj$direct, ncol = 1)
        Qobj$indirect <- matrix(Qobj$indirect, ncol = 1)
        Qobj$total <- matrix(Qobj$total, ncol = 1)
      }
      colnames(Qobj$direct) <- attr(x, "bnames")
      colnames(Qobj$indirect) <- attr(x, "bnames")
      colnames(Qobj$total) <- attr(x, "bnames")
      rownames(Qobj$direct) <- paste0("Q", 1:nrow(Qobj$direct))
      rownames(Qobj$indirect) <- paste0("Q", 1:nrow(Qobj$indirect))
      rownames(Qobj$total) <- paste0("Q", 1:nrow(Qobj$total))
      attr(mat, "Qobj") <- Qobj
    }
  }
  mat
}
