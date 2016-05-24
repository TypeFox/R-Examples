#' Report the total and relative number of each allele in a genotype object.
#'
#' @param genos1 A genotypes list object.
#'
#' @return Console output of the total an relative number of each allele.
#'
#' @examples \dontrun{reportGenos(preImputation)}

#' @export
reportGenos <- function(genos1) {

  cat(paste("The absolute number of genotypes in", deparse(substitute(genos1)),"is:"))
  t1 <- table(genos1$ABHmatrix)
  print(t1)

  cat(paste("\n"))
  cat(paste("or in percentage\n"))
  cat(paste("\n"))

  df1 <- data.frame("A" = round(100 / sum(t1) * t1[1], 3),
                    "B" = round(100 / sum(t1) * t1[2], 3),
                    "H" = round(100 / sum(t1) * t1[3], 3),
                    "N" = round(100 / sum(t1) * t1[4], 3))
  print(df1)

}
