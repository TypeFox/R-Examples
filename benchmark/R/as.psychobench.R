

#' Coerce benchmark experiment warehouse to preference table
#'
#' @param x A \code{\link{warehouse}} object
#' @param comparisons Return preference or performance table
#'
#' @return Data.frame with preference or performance table
#'
#' @export
as.psychobench <- function(x, comparisons = TRUE) {
  stopifnot(is(x, "warehouse"))


  ## Characteristics:
  ch <- x$viewDatasetCharacterization()
  stopifnot(nrow(ch) > 0)

  #ch <- subset(ch, samples != "basis")
  ch <- ch[ch$samples != "basis", ]
  ch$samples <- ch$samples[, drop = TRUE]

  ch <- reshape(ch, direction = "wide", v.names = "value",
                timevar = "characteristics",
                idvar = c("datasets", "samples"))
  attr(ch, "reshapeWide") <- NULL
  colnames(ch) <- sub("value.", "", colnames(ch))


  ## Performances:
  ap <- x$viewAlgorithmPerformance()
  stopifnot(nrow(ap) > 0)

  ap <- reshape(ap, direction = "wide", v.names = "value",
                timevar = c("algorithms"),
                idvar = c("datasets", "samples", "performances"))
  attr(ap, "reshapeWide") <- NULL
  colnames(ap) <- sub("value.", "", colnames(ap))


  if ( comparisons ) {
    ## Preference table:
    #pc <- subset.data.frame(ap, select = -c(samples, datasets, performances))
    pc <- ap[, -match(c("samples", "datasets", "performances"), names(ap))]
    pc <- bttree_paircomp(pc)

    #ret <- subset.data.frame(ch, select = -c(samples, datasets))
    ret <- ch[, -match(c("samples", "datasets"), names(ch))]
    ret$preference <- pc
  } else {
    ## Performance table:
    ret <- merge(ch, ap, sort = FALSE)
    #ret <- subset.data.frame(ret, select = -c(samples, datasets, performances))
    ret <- ret[, -match(c("samples", "datasets", "performances"), names(ret))]
  }

  class(ret) <- "data.frame"

  ret
}



#' @importFrom psychotools paircomp
bttree_paircomp <- function(x) {
  comprow <- function(x, tri) {
      eq <- outer(x, x, '==')[tri]
      g <- outer(x, x, '>')[tri]

      r <- as.numeric(!eq)
      r[!eq & g] <- 1
      r[!eq & !g] <- -1

      r
  }

  tri <- upper.tri(matrix(nrow = ncol(x),
                          ncol = ncol(x)))

  r <- t(apply(x, 1, comprow, tri))

  psychotools::paircomp(r, labels = colnames(x), mscale = c(-1, 0, 1))
}

