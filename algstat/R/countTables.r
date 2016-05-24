#' Count Similarly Margined Contingency Tables
#'
#' Count the number of contingency tables with the same marginals as a given table.
#'
#' \code{countTables} uses LattE's count function (via algstat's \code{\link{count}} function) to count the tables.  In many cases, the number of such tables is enormous.  In these cases, instead of giving back an integer \code{countTables} provides a character string with the integer in it; see examples.
#' 
#' @param table the table of interest
#' @param margins the margins to be fixed
#' @param dir directory to place the files in, without an ending /
#' @param opts options for count
#' @param quiet show latte output
#' @return an integer
#' @seealso \code{\link{count}}
#' @export countTables
#' @examples
#' \dontrun{ 
#' 
#' 
#' data(politics)
#' countTables(politics)
#'
#' data(handy)
#' countTables(handy)
#'
#' data(HairEyeColor)
#' eyeHairColor <- margin.table(HairEyeColor, 2:1)
#' countTables(eyeHairColor)
#' 
#' library(gmp)
#' as.bigz(countTables(eyeHairColor))
#' 
#'
#' # notice that even tables with small cells can have 
#' # huge fibers
#' data(drugs)
#' countTables(drugs)
#'
#' countTables(eyeHairColor, quiet = FALSE)
#' 
#'
#'
#' }
#' 
countTables <- function(table, 
    margins = as.list(1:length(dim(table))), dir = tempdir(), 
    opts = "", quiet = TRUE
){
  A <- hmat(dim(table), margins)
  cellNames <- paste0("t", colnames(A))
  
  margConds <- unname(apply(A, 1, function(v){
    paste(paste(v, cellNames), collapse = " + ")
  }))
  marginals <- as.integer(A %*% tab2vec(table))
  margConds <- paste0(margConds, " == ", marginals)
  
  nonnegConds <- paste0(cellNames, " >= 0")
  
  count(c(margConds, nonnegConds), dir, opts, quiet)  
}



