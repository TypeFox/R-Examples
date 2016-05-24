#' @title Row Percents

#' @description Computes row percentages for a given twoway table.
#' 
#' @rdname rowPerc
#' @usage rowPerc(tab)
#' @param tab A table, e.g.,
#' the result of \code{xtabs(~var1+var2,data=DataFrame)}.  
#' @return An object of class \code{table}, giving row percentages
#' for the input table.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' data(ledgejump)
#' MyTable <- xtabs(~weather+crowd.behavior,data=ledgejump)
#' rowPerc(MyTable)
rowPerc <-
function(tab)  {
  if (length(dim(tab))>1) {#tab is a two-way table
  rperc <- round(100*apply(tab,2,function(y)  y/rowSums(tab)),2)
  rperc <- cbind(rperc,rep(100,nrow(tab)))
  rownames(rperc) <- rownames(tab)
  colnames(rperc) <- c(colnames(tab),"Total")
  rperc2 <- as.table(rperc)
  names(dimnames(rperc2)) <- names(dimnames(tab))
  return(rperc2)
  } else {
    rperc <- round(100*tab/sum(tab),2)
    rperc <- as.matrix(rperc)
    rperc2 <- rbind(rperc,100)
    rperc2 <- t(rperc2)
    colnames(rperc2) <- c(rownames(rperc),"Total")
    rownames(rperc2) <- ""
    names(dimnames(rperc2)) <- c(names(dimnames(tab)),"")
    return(rperc2)
  }
}
