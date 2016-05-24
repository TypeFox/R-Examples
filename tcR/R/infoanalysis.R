#' Normalised log assymetry.
#' 
#' @description
#' Compute the value of the normalised log assymetry measure for the given data.frames
#' of the counts of shared clones.
#'
#' @param .alpha First mitcr data.frame or a list with data.frames.
#' @param .beta Second mitcr data.frame or NULL if \code{.alpha} is a list.
#' @param .by Which column use to merge. See "Details".
#' 
#' @details
#' Merge two data frames by the given column and compute
#' value Sum(Log((Percentage for shared clone S from alpha) / (Percentage for shared clone S from beta))) / (# of shared clones).
#' 
#' @return Value of the normalised log assymetry measure for the given data.frames.
assymetry<-function(.alpha, .beta = NULL, .by = 'CDR3.nucleotide.sequence'){
  if (class(.alpha) == 'list') {
    return(apply.asymm(.alpha, assymetry, .by = .by))
  }
  m<-merge(.alpha, .beta, by=.by)
  sum(log(m$Percentage.x/m$Percentage.y))/nrow(m)
}


#' Repertoires' analysis using information measures applied to V- and J- segment frequencies.
#' 
#' @aliases entropy.seg js.div.seg
#' 
#' @description
#' Information approach to repertoire analysis. Function \code{entropy.seg} applies Shannon entropy to V-usage and hence measures variability of V-usage.
#' Function \code{js.div.seg} applied Jensen-Shannon divergence to V-usage of two or more data frames and hence measures distance among this V-usages.
#' 
#' @usage
#' entropy.seg(.data, .genes = HUMAN_TRBV, .frame = c('all', 'in', 'out'),
#'             .quant = c(NA, "read.count", "umi.count", "read.prop", "umi.prop"),
#'             .ambig = F)
#' 
#' js.div.seg(.data, .genes = HUMAN_TRBV, .frame = c('all', 'in', 'out'),
#'            .quant = c(NA, "read.count", "umi.count", "read.prop", "umi.prop"),
#'            .norm.entropy = T, .ambig = F, .verbose = F, .data2 = NULL)
#' 
#' @param .data Mitcr data.frame or a list with mitcr data.frames.
#' @param .data2 NULL if .data is a list, or a second mitcr data.frame.
#' @param .genes Parameter to the \code{geneUsage} function.
#' @param .frame Character vector of length 1 specified which *-frames should be used:
#' only in-frame ('in'), out-of-frame ('out') or all sequences ('all').
#' @param .norm.entropy if T then divide result by mean entropy of 2 segments' frequencies. 
#' @param .ambig Parameter passed to \code{geneUsage}.
#' @param .quant Which column to use for the quantity of clonotypes: NA for computing only number of genes without using clonotype counts, "read.count" for the "Read.count" column, 
#' "umi.count" for the "Umi.count" column, "read.prop" for the "Read.proportion" column, "umi.prop" for 
#' the "Umi.proportion" column.
#' @param .verbose If T than output the data processing progress bar.
#' 
#' @return For \code{entropy.seg} - numeric integer with entropy value(s). For \code{js.div.seg} - integer of vector one if \code{.data} and \code{.data2} are provided;
#' esle matrix length(.data) X length(.data) if \code{.data} is a list.
#' 
#' @seealso \link{vis.heatmap}, \link{vis.group.boxplot}, \link{geneUsage}
entropy.seg <- function (.data, .genes = HUMAN_TRBV, .frame = c('all', 'in', 'out'),
                         .quant = c(NA, "read.count", "umi.count", "read.prop", "umi.prop"),
                         .ambig = F) {  
  if (class(.data) == 'list') {
    return(sapply(.data, entropy.seg, .quant = .quant, .frame = .frame, .genes = .genes, .ambig = .ambig))
  }
  
  .data <- get.frames(.data, .frame)
  
  if (has.class(.genes, "list") && length(.genes) == 2) {
    entropy(geneUsage(.data, .genes = .genes, .quant = .quant, .ambig = .ambig))
  } else {
    entropy(as.matrix(geneUsage(.data, .genes = .genes, .quant = .quant, .ambig = .ambig)[,-1]))
  }
}

js.div.seg <- function (.data, .genes = HUMAN_TRBV, .frame = c('all', 'in', 'out'),
                        .quant = c(NA, "read.count", "umi.count", "read.prop", "umi.prop"), .norm.entropy = T,
                        .ambig = F, .verbose = F, .data2 = NULL) {  
  if (class(.data) == 'list') {
    if (length(.data) == 2) {
      return(js.div.seg(.data[[1]], .genes, .frame, .quant, .norm.entropy, .ambig, .verbose, .data[[2]]))
    } else {
      return(apply.symm(.data, function (x, y) { js.div.seg(.data = x, .data2 = y, .quant = .quant, .frame = .frame, .ambig = .ambig, .norm.entropy = .norm.entropy, .genes = .genes) }, .verbose = .verbose))
    }
    
  }
  
  .data <- get.frames(.data, .frame)
  .data2 <- get.frames(.data2, .frame)
  
  if (has.class(.genes, "list") && length(.genes) == 2) {
    freq.alpha <- geneUsage(.data, .genes = .genes, .quant = .quant, .ambig = .ambig)
    freq.beta <- geneUsage(.data2, .genes = .genes, .quant = .quant, .ambig = .ambig)
  } else {
    freq.alpha <- geneUsage(.data, .genes = .genes, .quant = .quant, .ambig = .ambig)[,-1]
    freq.beta <- geneUsage(.data2, .genes = .genes, .quant = .quant, .ambig = .ambig)[,-1]
  }
  
  nrm = if (.norm.entropy) 0.5 * (entropy(freq.alpha) + entropy(freq.beta)) else 1
  js.div(freq.alpha, freq.beta) / nrm
}