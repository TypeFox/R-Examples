#' getCANSIM
#'
#' Extracts a complete CANSIM (Statistics Canada) data table
#' and converts it into a readily usable panel (wide) format.
#'
#' Geographic variables are renamed i, time variables are renamed t,
#' and all the other variables are renamed with a generic V1, V2, ..., Vn.
#' The generic variables keep the full Statistics Canada description by using a label.
#' @import reshape2 Hmisc utils
#'
#' @param cansimTableNumber - the table number we wish to retrieve from CANSIM.
#' @param showLabels - show the Statistics Canada labels after finishing extracting and converting the table, TRUE by default.
#' @return data frame containing CANSIM table.
#' @examples
#' getCANSIM(4010004, showLabels = FALSE)
#' getCANSIM(4010040)
#' getCANSIM(1530114)
#' @export
getCANSIM <- function(cansimTableNumber, showLabels = TRUE){
  df <- downloadCANSIM(cansimTableNumber)

  df2 <- dcast(df, df[,1] + df[,2] ~ StatCanVariable, value.var = "Value") #function from reshape2 package

  df2 <- df2[order(df2[,2]),]
  colnames(df2)[1] <- "t"
  colnames(df2)[2] <- "i"

  df3 <- labelCANSIM(df2)

  if(showLabels == TRUE) print( label(df3) )

  return(df3)
}
