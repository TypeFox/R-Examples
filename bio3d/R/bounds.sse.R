#' Obtain A SSE Object From An SSE Sequence Vector
#'
#' Inverse process of the funciton \code{\link{pdb2sse}}.
#'
#' @details call for its effects.
#'
#' @param x a character vector indicating SSE for each amino acid residue.
#' @param pdb an object of class \code{pdb} as obtained from
#'    function \code{\link{read.pdb}}. Can be ignored if \code{x} has 'names'
#'    attribute for residue labels.
#'
#' @return a 'sse' object. 
#'
#' @note In both \code{$helix} and \code{$sheet}, an additional 
#'    \code{$id} component is added to indicate the original numbering of the sse.
#'    This is particularly useful in e.g. \code{trim.pdb()} function.
#'
#' @seealso \code{\link{pdb2sse}}
#'
#' @author Xin-Qiu Yao & Barry Grant
#' 
#' @examples
#' \donttest{
#'    pdb <- read.pdb("1a7l")
#'    sse <- pdb2sse(pdb)
#'    sse.ind <- bounds.sse(sse)
#'    sse.ind 
#' } 
bounds.sse <- function(x, pdb=NULL) {

  if(length(x) == 0) return (NULL)
  
  strings <- names(x) 
  if(is.null(strings)) {
    if(!is.null(pdb)) {
       strings <- paste(pdb$atom[pdb$calpha, "resno"], pdb$atom[pdb$calpha, "chain"], 
                        pdb$atom[pdb$calpha, "insert"], sep = "_")
       if(length(strings) != length(x)) 
          stop("pdb doesn't match x")
    } else {
       strings <- paste(seq_along(x), NA, NA, sep="_")
    }
  } else {
    if(!is.null(pdb)) {
      warning("The x has 'names' attributes. The pdb is ignored")
    }
  }
   
  lstrings <- strsplit(strings, split="_")
  resno <- as.numeric(sapply(lstrings, "[", 1))
  chain <- sapply(lstrings, "[", 2)
  chain[chain == "NA"] <- ""
  insert <- sapply(lstrings, "[", 3)
  insert[insert == "NA"] <- ""
  id <- as.numeric(sapply(lstrings, "[", 4))

  sse.string <- paste(x, chain, id, sep="_")

# bounds doesn't work, use rle2 instead
  rl <- rle2(sse.string)
  inds <- cbind(seq_along(rl$inds), 
                start = c(1, rl$inds[-length(rl$inds)]+1), 
                end = rl$inds, length = rl$lengths )
#  inds <- bounds(sse.string, dup.inds=TRUE, pre.sort=FALSE) 

  # sort segments based on sse id (i.e. keep the original order of sse)
  ind.order <- order(id[inds[, "start"]])
  inds <- inds[ind.order, , drop = FALSE]

  # helix
  h.inds <- which(x[inds[, "start"]] == "H")
  if(length(h.inds) > 0) {

     h.id <- id[inds[h.inds, "start"]]
     if(any(is.na(h.id))) h.id <- seq_along(h.id)
 
     h <- list(start = resno[inds[h.inds, "start"]], end = resno[inds[h.inds, "end"]], 
        chain = chain[inds[h.inds, "start"]], id = h.id)
     names(h$start) <- insert[inds[h.inds, "start"]]
     names(h$end) <- insert[inds[h.inds, "end"]]
  } else {
     h <- list(start=NULL, end=NULL, chain=NULL, id=NULL)
  }

  # sheet
  e.inds <- which(x[inds[, "start"]] == "E")
  if(length(e.inds) > 0) {

     e.id <- id[inds[e.inds, "start"]]
     if(any(is.na(e.id))) e.id <- seq_along(e.id)
 
     e <- list(start = resno[inds[e.inds, "start"]], end = resno[inds[e.inds, "end"]], 
        chain = chain[inds[e.inds, "start"]], id = e.id)
     names(e$start) <- insert[inds[e.inds, "start"]]
     names(e$end) <- insert[inds[e.inds, "end"]]
  } else {
     e <- list(start=NULL, end=NULL, chain=NULL, id=NULL)
  }
 
  sse <- list(helix = h, sheet = e)
  class(sse) <- 'sse'  
  return( sse )
}
