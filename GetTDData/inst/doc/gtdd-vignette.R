## ----example1------------------------------------------------------------
library(GetTDData)

asset.code <- 'LTN'   # Identifier of assets 
maturity <- '010116'  # Maturity date as string (ddmmyy)

my.flag <- download.TD.data(asset.code = asset.code)

my.df <- read.TD.files(maturity = maturity)


## ----plot.prices, fig.width=7, fig.height=2.5----------------------------
library(ggplot2)

p <- ggplot(data = my.df, aes(x = as.Date(ref.date), y = price.bid, color = asset.code))
p <- p + geom_line(size = 1) + scale_x_date() + labs(title = '', x = 'Dates')
print(p)

## ----plot.yield, fig.width=7, fig.height=2.5-----------------------------
p <- ggplot(data = my.df, aes(x = as.Date(ref.date), y = yield.bid, color = asset.code))
p <- p + geom_line(size = 1) + scale_x_date()+ labs(title = '', x = 'Dates' )
print(p)


## ----example2, fig.width=7, fig.height=2.5-------------------------------
library(GetTDData)
library(ggplot2)

asset.code <- 'LTN'   # Name of asset
maturity <- NULL      # = NULL, downloads all maturities

# download data
my.flag <- download.TD.data(asset.code = asset.code)

# reads data
my.df <- read.TD.files(maturity = maturity)

# plot data (prices)
p <- ggplot(data = my.df, aes(x = as.Date(ref.date), y = price.bid, color = asset.code))
p <- p + geom_line() + scale_x_date() + labs(title = '', x = 'Dates', y = 'Prices' )
print(p)

# plot data (yields)
p <- ggplot(data = my.df, aes(x = as.Date(ref.date), y = yield.bid, color = asset.code))
p <- p + geom_line() + scale_x_date() + labs(title = '', x = 'Dates', y = 'Yields' )
print(p)


