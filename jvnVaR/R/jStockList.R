#library(jvnVaR)
jStockList <-
function(){
#data("stockList", envir=environment())
stockList <- NULL 
data('stockList', package='jvnVaR', envir=environment()) 
#data(stockList)
return(stockList)
}
