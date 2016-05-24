quoteStockXtsData <- function(x, ...){
    stock.df <- quoteStockTsData(x, ...)
    toXts(stock.df)
}
