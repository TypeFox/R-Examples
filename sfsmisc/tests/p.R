#### Plots etc
library(sfsmisc)

## A time-series with start and end *not* at year boundary:
data(EuStockMarkets)
SMI <- EuStockMarkets[, "SMI"]

p.ts(SMI)# gave warning (and was 'wrong' but "only" visually)
