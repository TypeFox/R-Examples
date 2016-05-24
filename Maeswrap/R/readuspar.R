#' Read understorey simulations
#' @description Reads the point-wise output file when the understorey was simulated.
#' @param filename The understorey file
#' @rdname plotuspar
#' @export
readuspar <- function(filename="uspar.dat"){

  watlines <- readLines(filename, 100)
  colloc <- grep("Columns",watlines)
  namesline <- watlines[colloc]
  NAMES <- delempty(strsplit(strsplit(namesline, ":")[[1]][2], " ")[[1]])
  uspar <- read.table(filename, header=FALSE, na.strings="-999.0000", skip=colloc)
  names(uspar) <- NAMES

return(uspar)
}
