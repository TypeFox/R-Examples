#' @name labelHaplotypes
#' @title Find and label haplotypes
#' @description Identify and group sequences that share the same haplotype.
#'
#' @param x a \code{\link{DNAbin}} \code{\link{multidna}}, or
#'   \linkS4class{gtypes} object.
#' @param prefix a character string giving prefix to be applied to numbered
#'   haplotypes. If NULL, haplotypes will be labeled with the first label
#'   from original sequences.
#' @param use.indels logical. Use indels when comparing sequences?
#' @param ... arguments to be passed to \code{labelHaplotypes.default}.
#'
#' @details If any sequences contain ambiguous bases (N's) they are first
#'   removed. Then haplotypes are assigned based on the remaining
#'   sequences. The sequences with N's that were removed are then assigned to
#'   the new haplotypes if it can be done unambiguously (they match only one
#'   haplotype with 0 differences once the N's have been removed). If this
#'   can't be done they are assigned NAs and listed in the
#'   \code{unassigned} element.
#'
#' @return
#'   \code{DNAbin} or \code{multidna}, a list with the following elements:
#'     \tabular{ll}{
#'       \code{haps} \tab named vector (\code{DNAbin}) or list of named vectors
#'         (\code{multidina}) of haplotypes for each sequence in \code{x}.\cr
#'       \code{hap.seqs} \tab \code{DNAbin} or \code{multidna} object
#'         containing sequences for each haplotype.\cr
#'       \code{unassigned} \tab \code{data.frame} listing closest matching
#'         haplotypes and the number of substitutions different. Will be
#'         \code{NULL} if no sequences remain unassigned.
#'     }
#'   \code{gtypes}, a list with the following elements:
#'     \tabular{ll}{
#'       \code{gtypes} \tab the new \code{gtypes} object with the haplotypes
#'         reassigned.\cr
#'       \code{unassigned} \tab a list containing the \code{unassigned}
#'         attribute \code{data.frame} for each gene if present,
#'         otherwise \code{NULL}.\cr
#'      }
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' # create 5 example short haplotypes
#' haps <- c(
#'   H1 = "ggctagct",
#'   H2 = "agttagct",
#'   H3 = "agctggct",
#'   H4 = "agctggct",
#'   H5 = "ggttagct"
#' )
#
#' # draw and label 100 samples
#' sample.seqs <- sample(names(haps), 100, rep = TRUE)
#' ids <- paste(sample.seqs, 1:length(sample.seqs), sep = "_")
#' sample.seqs <- lapply(sample.seqs, function(x) strsplit(haps[x], "")[[1]])
#' names(sample.seqs) <- ids
#'
#' # add some random ambiguities
#' with.error <- sample(1:length(sample.seqs), 10)
#' for(i in with.error) {
#'   site <- sample(1:length(sample.seqs[[i]]), 1)
#'   sample.seqs[[i]][site] <- "n"
#' }
#'
#' # convert to DNAbin
#' library(ape)
#' sample.seqs <- as.DNAbin(sample.seqs)
#'
#' hap.assign <- labelHaplotypes(sample.seqs, prefix = "Hap.")
#' hap.assign
#'
#' @importFrom swfscMisc zero.pad
#' @export
#'
labelHaplotypes <- function(x, prefix = NULL, use.indels = FALSE) {
  UseMethod("labelHaplotypes")
}

#' @rdname labelHaplotypes
#' @export
#'
labelHaplotypes.default  <- function(x, prefix = NULL, use.indels = TRUE) {
  if(!inherits(x, "DNAbin")) stop("'x' must be a DNAbin object.")
  x <- as.matrix(x)
  
  # return same data if only one sequence exists
  if(nrow(x) == 1) {
    haps <- rownames(x)
    names(haps) <- haps
    return(list(haps = haps, hap.seqs = x, unassigned = NULL))
  }
  
  # find sequences without Ns
  has.ns <- apply(as.character(x), 1, function(bases) "n" %in% tolower(bases))
  if(sum(!has.ns) == 1) {  
    warning("There is only one sequence without ambiguities (N's). Can't assign haplotypes. NULL returned.",
            call. = FALSE, immediate. = TRUE)
    return(NULL)
  }

  # get pairwise distances and set all non-0 distances to 1
  x.no.ns <- x[!has.ns, ]
  hap.dist <- dist.dna(x.no.ns, model = "N", pairwise.deletion = TRUE)
  if(use.indels) hap.dist <- hap.dist + dist.dna(x.no.ns, model = "indelblock")
  hap.dist <- as.matrix(hap.dist)
  hap.dist[hap.dist > 0] <- 1

  # create haplotype code out of 0s and 1s
  hap.code <- as.numeric(factor(apply(hap.dist, 1, paste, collapse = "")))
  names(hap.code) <- rownames(hap.dist)

  # rename haplotypes
  hap.labels <- if(!is.null(prefix)) {
    # use prefix+number if prefix given
    # sort based on frequency first
    hap.order <- as.numeric(names(sort(table(hap.code), decreasing = TRUE)))
    hap.nums <- zero.pad(1:length(hap.order))
    names(hap.order) <- paste(prefix, hap.nums, sep = "")
    names(sort(hap.order))
  } else {
    # if no prefix, use first sequence name for each haplotype
    #hap.code.sort <- hap.code[order(names(hap.code))]
    #names(sort(hap.code[!duplicated(hap.code.sort)]))
    names(hap.code[!duplicated(hap.code)])
  }
  hap.code <- hap.labels[hap.code]
  names(hap.code) <- rownames(hap.dist)

  # get sequences for each haplotype
  unique.codes <- hap.code[!duplicated(hap.code)]
  hap.seqs <- x[names(unique.codes), , drop = FALSE]
  rownames(hap.seqs) <- unique.codes
  hap.seqs <- hap.seqs[order(rownames(hap.seqs)), , drop = FALSE]

  unassigned.df <- if(!all(!has.ns)) {
    # calculate distance between haplotypes and samples with Ns
    with.ns.seqs <- x[has.ns, ]
    haps.and.ns <- rbind(hap.seqs, with.ns.seqs)
    hap.dist <- dist.dna(haps.and.ns, model = "N", pairwise.deletion = TRUE)
    if(use.indels) hap.dist <- hap.dist + dist.dna(haps.and.ns, model = "indelblock")
    hap.dist <- as.matrix(hap.dist)[rownames(with.ns.seqs), rownames(hap.seqs), drop = FALSE]

    # match samples with Ns to one or more haplotypes
    ns.matches <- apply(hap.dist, 1, function(d) {
      min.dist <- min(d)
      haps <- names(d)[d == min.dist]
      list(haps = haps, min.dist = min.dist)
    })

    is.matched <- sapply(ns.matches, function(match.info) {
      match.info$min.dist == 0 & length(match.info$haps) == 1
    })

    # assign matched haplotypes
    assigned <- lapply(ns.matches[is.matched],
                       function(match.info) match.info$haps)
    assigned <- do.call(c, assigned)

    # identify haplotypes that can't be assigned
    unassigned.df <- lapply(ns.matches[!is.matched], function(match.info) {
      data.frame(closest.match = paste(match.info$haps, collapse = ", "),
                 dist = match.info$min.dist)
    })
    unassigned.df <- do.call(rbind, unassigned.df)
    unassigned <- if(is.null(unassigned.df)) {
      NULL
    } else {
      warning("Some sequences could not be unambiguously assigned to a haplotype:",
              call. = FALSE, immediate. = TRUE)
      print(unassigned.df)
      cat("\n")
      rep(NA, nrow(unassigned.df))
    }
    names(unassigned) <- rownames(unassigned.df)

    # create full haplotype assignment vector
    hap.code <- c(hap.code, assigned, unassigned)[rownames(x)]
    unassigned.df
  } else NULL

  list(haps = hap.code, hap.seqs = hap.seqs, unassigned = unassigned.df)
}

#' @rdname labelHaplotypes
#' @export
#'
labelHaplotypes.gtypes <- function(x, ...) {
  # check that sequences are present
  if(ploidy(x) > 1 | is.null(sequences(x))) {
    stop("'x' is not haploid or does not have any sequences")
  }
  
  # label haplotypes for each gene
  new.haps <- lapply(
    getSequences(sequences(x), simplify = FALSE), labelHaplotypes, ...
  )
  has.errors <- sapply(new.haps, is.null)
  if(sum(has.errors) > 0) {
    msg <- paste(names(new.haps)[has.errors], collapse = ", ")
    msg <- paste("haplotypes could not be assigned for:", msg)
    stop(msg)
  }
  
  # create haplotype data.frame
  hap.df <- gtypes2df(x)
  for(gene in names(new.haps)) {
    old.haps <- hap.df[, gene]
    hap.df[, gene] <- new.haps[[gene]]$haps[old.haps]
  }
  
  # collect sequences
  hap.seqs <- lapply(new.haps, function(x) x$hap.seqs)
  names(hap.seqs) <- names(new.haps)
  
  # collect unassigned
  unassigned <- lapply(new.haps, function(x) x$unassigned)
  
  # create new gtypes
  st <- strata(x)
  x <- df2gtypes(
    hap.df, ploidy = 1, id.col = 1, strata.col = 2, loc.col = 3,
    sequences = hap.seqs, description = description(x), schemes = schemes(x), 
    other = other(x)
  )
  strata(x) <- st

  list(gtypes = x, unassigned = unassigned)
}
