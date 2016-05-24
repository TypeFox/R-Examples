#' General function for the repertoire diversity estimation.
#' 
#' @description
#' General interface to all cloneset diversity functions.
#' 
#' @param .data Cloneset or a list of clonesets.
#' @param .method Which method to use for the diversity estimation. See "Details" for methods.
#' @param .quant Which column to use for the quantity of clonotypes: "read.count" for the "Read.count" column, 
#' "umi.count" for the "Umi.count" column, "read.prop" for the "Read.proportion" column, "umi.prop" for 
#' the "Umi.proportion" column.
#' @param .q q-parameter for the Diversity index.
#' @param .norm If T than compute the normsalised entropy.
#' @param .do.norm One of the three values - NA, T or F. If NA than check for distrubution (sum(.data) == 1)
#' and normalise it with the given laplace correction value if needed. if T then do normalisation and laplace
#' correction. If F than don't do normalisaton and laplace correction.
#' @param .laplace Value for Laplace correction.
#' 
#' @details
#' You can see a more detailed description for each diversity method at \link{diversity}.
#' 
#' Parameter \code{.method} can have one of the following value each corresponding to the specific method:
#' 
#' - "div" for the true diversity, or the effective number of types (basic function \code{diversity}).
#' 
#' - "inv.simp" for the inverse Simpson index (basic function \code{inverse.simpson}).
#' 
#' - "gini" for the Gini coefficient (basic function \code{gini}).
#' 
#' - "gini.simp" for the Gini-Simpson index (basic function \code{gini.simpson}).
#' 
#' - "chao1" for the Chao1 estimator (basic function \code{chao1}).
#' 
#' - "entropy" for the Shannon entropy measure (basic function \code{entropy}).
#' 
#' @seealso \link{diversity}, \link{entropy}
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' twb.div <- repDiversity(twb, "chao1", "read.count")
#' }
repDiversity <- function (.data,
                          .method = c("chao1", "gini.simp", "inv.simp", "gini", "div", "entropy"), 
                          .quant = c("read.count", "umi.count", "read.prop", "umi.prop"), 
                          .q = 5,
                          .norm = F,
                          .do.norm = NA,
                          .laplace = 0) {
  
  quant <- .column.choice(.quant, T)
  
  fun <- switch(.method[1], 
                chao1 = function (x, ...) chao1(x),
                gini.simp = gini.simpson,
                inv.simp = inverse.simpson,
                gini = gini,
                div = function (x, ...) diversity(x, .q = .q, ...),
                entropy = function (x, ...) entropy(x, .norm = .norm, ...),
                { .verbose.msg("You have specified an invalid method identifier. Choosed method: chao1\n", T); chao1 })
  
  if (has.class(.data, 'data.frame')) { .data <- list(Sample = .data) }
  
  .data <- .fix.listnames(.data)
  
  sapply(.data, function (x) fun(x[[quant]], .do.norm = .do.norm, .laplace = .laplace))
}