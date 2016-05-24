toXts <- function(stock.df){
    stock.matrix <- as.matrix(stock.df[,-1]);
    colnames(stock.matrix) <- c("Open", "High", "Low", "Close", "Volume", "AdjClose")
    rownames(stock.matrix) <- as.character(stock.df$date)
    stock.xts <- as.xts(stock.matrix, descr='stock xts data')

    stock.xts
}