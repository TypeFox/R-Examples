## ----message=F, warning=F, results='hide', echo=F------------------------
library(ggplot2)
library(scales)
library(bdscale)

## ------------------------------------------------------------------------
data(nyse)

## ------------------------------------------------------------------------
set.seed(12345)
df <- data.frame(date=nyse, price=cumsum(rnorm(length(nyse))) + 100)
df <- subset(df, as.Date('2014-08-01') <= date & date <= as.Date('2014-10-08'))

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(tail(df, 5))

## ------------------------------------------------------------------------
plot <- ggplot(df, aes(x=date, y=price)) + geom_step() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

## ---- fig.width=7--------------------------------------------------------
plot + ggtitle('calendar dates')

## ---- fig.width=7--------------------------------------------------------
plot + scale_x_bd(business.dates=nyse, labels=date_format("%b '%y")) + 
  ggtitle('business dates, month breaks')

## ---- fig.width=7--------------------------------------------------------
plot + scale_x_bd(business.dates=nyse, max.major.breaks=10, labels=date_format('%d %b')) + 
  ggtitle('business dates, week breaks')

## ---- fig.width=7--------------------------------------------------------
options <- as.Date(c('2014-08-15', '2014-09-19'))

plot + 
  geom_vline(xintercept=as.numeric(options), size=2, alpha=0.25) + 
  ggtitle('calendar dates, option expiry')

## ---- fig.width=7--------------------------------------------------------
plot + 
  geom_vline(xintercept=bd2t(options, business.dates=nyse), size=2, alpha=0.25) + 
  scale_x_bd(business.dates=nyse) +
  ggtitle('business dates, option expiry')

