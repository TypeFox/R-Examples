#' General function for the repertoire overlap evaluation.
#' 
#' @description
#' General interface to all cloneset overlap functions.
#' 
#' @param .data List of clonesets.
#' @param .method Which method to use for the overlap evaluation. See "Details" for methods.
#' @param .seq Which clonotype sequences to use for the overlap: "nuc" for "CDR3.nucleotide.sequence", "aa" for 
#' "CDR3.amino.acid.sequence".
#' @param .quant Which column to use for the quantity of clonotypes: "read.count" for the "Read.count" column, 
#' "umi.count" for the "Umi.count" column, "read.prop" for the "Read.proportion" column, "umi.prop" for 
#' the "Umi.proportion" column. Used in "morisita" and "horn".
#' @param .vgene If T than use V genes in computing shared or similar clonotypes. Used in all methods.
#' @param .norm If T than compute the normalised number of shared clonotypes. Used in "exact".
#' @param .a,.b Alpha and beta parameters for "tversky". Default values gives the Jaccard index measure.
#' @param .do.unique If T than remove duplicates from the input data, but add their quantities to their clones.
#' @param .verbose If T than output the data processing progress bar.
#' 
#' @details
#' You can see a more detailed description for each overlap method at \link{intersectClonesets} and \link{similarity}.
#' 
#' Parameter \code{.method} can have one of the following value each corresponding to the specific method:
#' 
#' - "exact" for the shared number of clonotypes (basic function \code{intersectClonesets(..., .type = "..e")}).
#' 
#' - "hamm" for the number of similar clonotypes by the Hamming distance (basic function \code{intersectClonesets(..., .type = "..h")}).
#' 
#' - "lev" for the number of similar clonotypes by the Levenshtein distance (basic function \code{intersectClonesets(..., .type = "..l")}).
#' 
#' - "jaccard" for the Jaccard index (basic function \code{jaccard.index}).
#' 
#' - "morisita" for the Morisita's overlap index (basic function \code{morisita.index}).
#' 
#' - "tversky" for the Tversky index (basic function \code{tversky.index}).
#' 
#' - "overlap" for the overlap coefficient (basic function \code{overlap.coef}).
#' 
#' - "horn" for the Horn's index (basic function \code{horn.index}).
#' 
#' @seealso  \link{intersectClonesets}, \link{similarity}, \link{repDiversity}
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' repOverlap(twb, "exact", .seq = "nuc", .vgene = F)
#' repOverlap(twb, "morisita", .seq = "aa", .vgene = T, .quant = "umi.count")
#' ov <- repOverlap(twb)
#' vis.pca(prcomp(ov, scale. = T), list(A = c(1, 2), B = c(3, 4)))
#' }
repOverlap <- function (.data,
                        .method = c("exact", "hamm", "lev", "jaccard", "morisita", "tversky", "overlap", "horn"), 
                        .seq = c("nuc", "aa"), 
                        .quant = c("read.count", "umi.count", "read.prop", "umi.prop"),
                        .vgene = F,
                        .norm = T,
                        .a = .5,
                        .b = .5,
                        .do.unique = T,
                        .verbose = T) {
  
  .merge.with.v <- function (.data, .seqcol, .vgene) {
    if (.vgene) {
      lapply(.data, function (x) paste0(x[[.seqcol]], x$V.gene) )
    } else {
      lapply(.data, function (x) x[[.seqcol]] )
    }
  }
  
  .pair.fun <- function (.data, .fun, .verbose) {
    if (length(.data) == 2) { .fun(.data[[1]], .data[[2]]) }
    else { apply.symm(.data, .fun, .verbose = .verbose) }
  }
  
  
  quant <- .column.choice(.quant, .verbose)
  
  if (!has.class(.data, 'list')) { cat("Error! Input data MUST be a list of clonesets!\n"); return(NA); }
  
  if (length(.data) < 2) { cat("Error! Input data MUST be a list of clonesets of minimum length 2!\n"); return(NA); }
  
  .data <- .fix.listnames(.data)
  
  .seqcol <- "CDR3.nucleotide.sequence"
  if (.seq[1] == "aa") { .seqcol <- "CDR3.amino.acid.sequence" }
  
  if (.method[1] %in% c("exact", "hamm", "lev")) {
    let1 <- "n"
    if (.seq[1] == "aa") { let1 <- "a" }
    let2 <- "0"
    if (.vgene) { let2 <- "v" }
    let3 <- substr(.method[1], 1, 1)
    .type <- paste0(let1, let2, let3)
    
    if (length(.data) != 2) {
      intersectClonesets(.data, .type = .type, .norm = .norm, .verbose = .verbose)
    } else {
      params <- list(.type = .type, .norm = .norm, .verbose = .verbose)
      do.call(intersectClonesets, c(list(.data[[1]], .data[[2]]), params))
    }
  }
  else if (.method[1] %in% c("morisita", "horn")) {
    .fun <- function (x, y) { horn.index(x, y, F) }
    if (.method[1] == "morisita") {
      .fun <- function (x, y) { morisitas.index(x, y, F) }
    }
    
    .verbose.msg("Preprocessing data...\t", .verbose)
    new.data <- .merge.with.v(.data, .seqcol, .vgene)
    new.reads <- lapply(.data, "[[", quant)
    
    if (.do.unique) {
      if (.verbose) { pb <- set.pb(length(.data)) }
      for (i in 1:length(.data)) {
        .data[[i]] <- as.data.frame(summarise(grouped_df(data.frame(Sequence = new.data[[i]], Count = new.reads[[i]], stringsAsFactors = F),
                                           list(as.name("Sequence"))), Count = sum(Count)), stringsAsFactors = F)
        if (.verbose) { add.pb(pb) }
      }
      if (.verbose) { close(pb) }
    }
    
    .pair.fun(.data,
              function (x, y) { .fun(x, y) },
              .verbose)
  }
  else if (.method[1] == "jaccard") {
    .pair.fun(.merge.with.v(.data, .seqcol, .vgene),
              function (x, y) { jaccard.index(x, y) },
              .verbose)
  }
  else if (.method[1] == "overlap") {
    .pair.fun(.merge.with.v(.data, .seqcol, .vgene),
              function (x, y) { overlap.coef(x, y) },
              .verbose)
  }
  else if (.method[1] == "tversky") {
    .pair.fun(.merge.with.v(.data, .seqcol, .vgene),
              function (x, y) { tversky.index(x, y, .a = .a, .b = .b) },
              .verbose)
  }
  else {
    .verbose.msg("You have specified an invalid method identifier. Please check your input arguments.\n", .verbose)
    return(NA)
  } 
}