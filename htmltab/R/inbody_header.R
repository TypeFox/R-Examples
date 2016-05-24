#' Reshape in table header information into wide format
#' @param tab the table data frame
#' @param table.Node the table node
#' @param trindex the tr index of the inbody rows
#' @param xpath the xpath for the inbody rows
#' @return the modified R data frame
#' @noRd
create_inbody <- function(tab, table.Node, trindex, xpath){

  # When no inbody header was specified
  if(length(trindex) == 0){
    return(tab)
  }

  # Else ...
  N <- length(trindex)
  N.row <- nrow(tab)
  df <- data.frame(matrix(NA, ncol = N, nrow = N.row))
  colnames(df) <- paste0("Header_", 1:N)

  for(i in 1:N){
    val <- XML::xpathSApply(table.Node, xpath[[i]], XML::xmlValue)
    index <- trindex[[i]]
    index <- c(index, (N.row + 1))

    df[index[1]:(N.row), (N+1-i)] <- rep(val, diff(index))
  }

  # Combine
  tab <- cbind(df, tab)

  return(tab)
}
