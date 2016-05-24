#' @title Column Percents

#' @description Computes column percentages for a given twoway table.
#' 
#' @rdname colPerc
#' @usage colPerc(tab)
#' @param tab A two way table, e.g.,
#' the result of \code{xtabs(~var1+var2,data=DataFrame)}.  
#' @return An object of class \code{table}, giving column percentages
#' for the input table.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' MyTable <- xtabs(~weather+crowd.behavior,data=ledgejump)
#' colPerc(MyTable)
colPerc <- function (tab) 
{
  cperc <- t(round(100 * apply(tab, 1, function(y) y/colSums(tab)), 
                   2))
  cperc <- rbind(cperc, rep(100, ncol(tab)))
  colnames(cperc) <- colnames(tab)
  rownames(cperc) <- c(rownames(tab), "Total")
  names(dimnames(cperc)) <- names(dimnames(tab))
  cperc
}
