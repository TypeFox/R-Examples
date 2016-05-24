keepUniqWell.fn <-
function( dataIn.df, colFocus=c("Barcode","Rowpos","Colpos") )
{
#*****************************************************************#
# author: Xiaohua Douglas Zhang, 2006                             #     
# keep rows with unique wells in each plate, without side effects #  
#*****************************************************************#
  uniqData = unique(dataIn.df[, colFocus])
  uniqNames = dimnames(uniqData)[[1]]
  dataOut.df = dataIn.df[uniqNames,]
  dataOut.df
}
