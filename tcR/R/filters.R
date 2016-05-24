########## Various ways to subset the data ##########


#' Contamination filtering.
#' 
#' @aliases contamination.stats decontamination
#' 
#' @description 
#' Occasionally DNA or RNA libraries are contaminate each other. To address this issue and estimate contamination rate \code{tcR} offers 
#' \code{contamination.stats} and \code{decontamination} functions. The \code{decontamination} function received data
#'  (either data frame or a list with data frames) and a limit for clonal proportion as arguments. 
#'  Script searches for a similar clones to the first data frame in the other (or performs pairwise searches if the given data is a list) 
#'  and removes clones from the first data frame, which has been found in the second one with counts less or equal to 10 * counts of similar clones 
#'  in the first one. Function \code{contamination.stats} will return the number of clones which will be removed with the \code{contamination.stats} function.
#' 
#' @usage
#' contamination.stats(.data1, .data2, .limit = 20, .col = 'Read.count')
#' 
#' decontamination(.data1, .data2, .limit = 20, .col = 'Read.count', .symm = T)
#' 
#' @param .data1 First data frame with columns 'CDR3.nucleotide.sequence' and 'Read.count'. Will be checked for contamination.
#' @param .data2 Second data frame with such columns. Will be used for checking for sequences which contaminated the first one.
#' @param .limit Parameter for filtering: all sequences from \code{.data1} which are presented in \code{.data2} and (count of  in \code{.data2}) / (count of seq in \code{.data1}) >= \code{.limit} are removed.
#' @param .col Column's name with clonal count.
#' @param .symm if T then perform filtering out of sequences in .data1, and then from .data2. Else only from .data1.
#' 
#' @return Filtered \code{.data1} or a list with filtered both \code{.data1} and \code{.data2}.
contamination.stats <- function (.data1, .data2, .limit = 20, .col = 'Read.count') {
  if (has.class(.data1, 'list') && has.class(.data2, 'list')) {
    res1 <- sapply(1:length(.data1), function (i) contamination.stats(.data1[[i]], .data2[[i]], .limit, .col))
    colnames(res1) <- paste0(names(.data1), ':', names(.data2))
    res2 <- sapply(1:length(.data1), function (i) contamination.stats(.data2[[i]], .data1[[i]], .limit, .col))
    colnames(res2) <- paste0(names(.data2), ':', names(.data1))
    return(cbind(res1, res2))
  }
  
  inds <- intersectIndices(.data1$CDR3.nucleotide.sequence, .data2$CDR3.nucleotide.sequence, 'exact')
  logic <- .data2[inds[,2], .col] / .data1[inds[,1], .col] >= .limit
  counts <- .data1[logic, .col]
  c(Nclones.del = length(counts), summary(counts))
}

decontamination <- function (.data1, .data2, .limit = 20, .col = 'Read.count', .symm = T) {
  .filter <- function (.data1, .data2, .limit, .col) {    
    inds <- intersectIndices(.data1$CDR3.nucleotide.sequence, .data2$CDR3.nucleotide.sequence, 'exact')
    logic <- .data2[inds[,2], .col] / .data1[inds[,1], .col] >= .limit
    .data1[!logic, ]
  }
  
  if (has.class(.data1, 'list') && has.class(.data2, 'list')) {
    res <- lapply(1:length(.data1), function (i) decontamination(.data1[[i]], .data2[[i]], .limit, .col, T))
    
    return(res)
  }
  
  if (.symm) {
    list(First = .filter(.data1, .data2, .limit, .col), Second = .filter(.data2, .data1, .limit, .col))
  } else {
    .filter(.data1, .data2, .limit, .col)
  }
}