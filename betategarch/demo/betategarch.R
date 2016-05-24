#load data, define zoo object, plot data:
data(nasdaq)
y <- zoo(nasdaq[,"nasdaqret"], order.by=as.Date(nasdaq[,"day"], "%Y-%m-%d"))
plot(y, main="The Nasdaq 100 index (daily)", xlab="", ylab="Log-return in %")

#estimate 1-comp model, plot standard deviations:
nasdaq1comp <- tegarch(y)
nasdaq1stdev <- fitted(nasdaq1comp)
plot(nasdaq1stdev, main="", ylab="1-comp: St.dev.", xlab="")
nasdaq1comp

#estimate 2-comp model, plot standard deviations:
nasdaq2comp <- tegarch(y, components=2)
nasdaq2stdev <- fitted(nasdaq2comp)
plot(nasdaq2stdev, main="", ylab="2-comp: St.dev.", xlab="")
nasdaq2comp

#predict standard deviations up to 5-steps ahead:
set.seed(123)
predict(nasdaq1comp, n.ahead=5)