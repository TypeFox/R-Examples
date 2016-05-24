##' @title Get Year Measurement
##' @description using column_name_measurement_date column in the form MM/YYYY creates a new column with the name "ANO_MEDICAO" in YYYY format
##' @param dataFrame that has the column DATE(MM/YYYY) and a ID column_name_plot
##' @param column_name_measurement_date column with a date format
##' @param column_name_plot a column of dataFrame, identification of plot (ID_plot)
##' @return dataFrame dataframe that has columns column_name_measurement_date, column_name_plot, ANO_MEDICAO
##' @import stringr
##' @import plyr
##' @import tcltk
##' @examples
##' column_name_measurement_date <- c("02/2009","02/2010","02/2011","02/2012")
##' column_name_plot <- c(1,2,3,4)
##' test <- data.frame(column_name_measurement_date,column_name_plot)
##' getAnoMedicao(test,"column_name_measurement_date","column_name_plot")
##' @export
getAnoMedicao <- function (dataFrame, column_name_measurement_date, column_name_plot) {

  #Get the name of dataFrame.
  nomedataFrame = toString(substitute(dataFrame))
  sql ="SELECT distinct CD_PARCELA, DATA_MEDICAO from dataFrame group by CD_PARCELA, DATA_MEDICAO"
  sql = gsub("dataFrame", nomedataFrame, sql)
  sql = gsub ("CD_PARCELA",column_name_plot,sql)
  sql = gsub ("DATA_MEDICAO",column_name_measurement_date,sql)

  df = NULL;
  dataFrame = sqldf(sql)
  comando = NULL;
  comando = "df <- ldply(str_split(dataFrame$DATA_MEDICAO,  '/'))[,2]"
  comando = gsub("DATA_MEDICAO",column_name_measurement_date,comando)
 eval(parse(text=(comando)))

  dataFrame = cbind(dataFrame, df)
  #id_col = which(names(dataFrame) == "df") # Get the identifier column of a dataframe.
  names(dataFrame)[3] = "ANO_MEDICAO"
  dataFrame$ANO_MEDICAO = as.numeric(str_c(dataFrame$ANO_MEDICAO))

  #to ensure it will not be concatenated with any date that already has 0 20 before.
  if(str_length(dataFrame$ANO_MEDICAO[1])==2)
  {
    dataFrame$ANO_MEDICAO = as.numeric(str_c("20", dataFrame$ANO_MEDICAO))
  }

  return (dataFrame)
}
