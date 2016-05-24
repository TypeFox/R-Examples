############################################################
# Part 1 - Define trading signals and construct portfolios
############################################################

require(PortfolioEffectHFT)

# Moving average
MA=function(x,order){
  result = x
  x1 = c(0,x)
  result[(order):NROW(x)] = (cumsum(x1)[-(1:(order))]-cumsum(x1)[-((NROW(x1)-order+1):NROW(x1))])/order
  result[1:(order-1)] = cumsum(x[1:(order-1)])/(1:(order-1))
  return(result-0.0000000001)
}

# Lag function padded with zeroes
lagpad=function(x, k=1) {
  i = is.vector(x)
  if(is.vector(x)) x = matrix(x) else x = matrix(x,nrow(x))
  if(k>0) {
    x = rbind(matrix(rep(0, k*ncol(x)),ncol=ncol(x)), matrix(x[1:(nrow(x)-k),], ncol=ncol(x)))
  }
  else {
    x = rbind(matrix(x[(-k+1):(nrow(x)),], ncol=ncol(x)),matrix(rep(0, -k*ncol(x)),ncol=ncol(x)))
  }
  if(i) x[1:length(x)] else x
}


symbol = "GOOG"
dateStart = "2014-10-13 09:30:00"
dateEnd = "2014-10-14 16:00:00"
portfolio = portfolio_create(dateStart,dateEnd)

portfolio_addPosition(portfolio,symbol,1)
price = position_price(portfolio,symbol)
printTime = price[,1]

# Create two strategy by moving average differents length
highFrequencyStrategy = array(0,dim=NROW(price))
highFrequencyStrategy[price[,"value"]<MA(price[,"value"],200)]=100
highFrequencyStrategy[(NROW(price)/2):NROW(price)]=0
highFrequencyStrategy=lagpad(highFrequencyStrategy)

lowFrequencyStrategy=array(0,dim=NROW(price))
lowFrequencyStrategy[price[,"value"]<MA(price[,"value"],800)]=100
lowFrequencyStrategy=lagpad(lowFrequencyStrategy)

############################################################
# Part 2 - Holding Times
############################################################

highFrequencyStrategyPlot = NULL
highFrequencyStrategyPlot = rbind(highFrequencyStrategyPlot,data.frame(positions=sum(highFrequencyStrategy>0)/60,Legends="Has position"))
highFrequencyStrategyPlot = rbind(highFrequencyStrategyPlot,data.frame(positions=sum(highFrequencyStrategy==0)/60,Legends="No position"))
lowFrequencyStrategyPlot = NULL
lowFrequencyStrategyPlot = rbind(lowFrequencyStrategyPlot,data.frame(positions=sum(lowFrequencyStrategy>0)/60,Legends="Has position"))
lowFrequencyStrategyPlot = rbind(lowFrequencyStrategyPlot,data.frame(positions=sum(lowFrequencyStrategy==0)/60,Legends="No position"))

xlabel = ""
ylabel = "In minutes"
p1 = ggplot(highFrequencyStrategyPlot, aes(x = "", y = positions, fill = Legends))+xlab(xlabel)+ylab(ylabel) +  geom_bar(stat = "identity")+
  geom_text(aes(x= 1,y = positions/2 + c(0, cumsum(positions)[-length(positions)]), label = paste(round(positions,digits =1)," minutes",sep="")), size=7,col="#d5e4eb")
p1 = p1+coord_polar("y")+ ggtitle("Intraday holding period for\n high-frequency strategy")+
  util_plotTheme()+util_fillScheme()
p2 = ggplot(lowFrequencyStrategyPlot, aes(x = "", y = positions, fill = Legends))+xlab(xlabel)+ylab(ylabel) +  geom_bar(stat = "identity")+
  geom_text(aes(x= 1,y = positions/2 + c(0, cumsum(positions)[-length(positions)]), label = paste(round(positions,digits =1)," minutes",sep="")), size=7,col="#d5e4eb")
p2 = p2+coord_polar("y")+ ggtitle("Intraday holding period for\n low-frequency strategy")+
  util_plotTheme()+util_fillScheme()
util_multiplot(p1,p2,cols=2)

############################################################
# Part 3 - Holding Intervals
############################################################

highFrequencyPortfolioHoldOnly = portfolio_create(dateStart,dateEnd)
portfolio_settings(highFrequencyPortfolioHoldOnly,holdingPeriodsOnly=TRUE)
portfolio_addPosition(highFrequencyPortfolioHoldOnly,symbol,quantity=as.numeric(highFrequencyStrategy),time=printTime)
highFrequencyPortfolioHoldOnly

highFrequencyPortfolioAllDay = portfolio_create(dateStart,dateEnd)
portfolio_settings(highFrequencyPortfolioAllDay,holdingPeriodsOnly=FALSE)
portfolio_addPosition(highFrequencyPortfolioAllDay,symbol,quantity=as.numeric(highFrequencyStrategy),time=printTime)
highFrequencyPortfolioAllDay

lowFrequencyPortfolioHoldOnly = portfolio_create(dateStart,dateEnd)
portfolio_settings(lowFrequencyPortfolioHoldOnly,holdingPeriodsOnly=TRUE)
portfolio_addPosition(lowFrequencyPortfolioHoldOnly,symbol,quantity=as.numeric(lowFrequencyStrategy),time=printTime)
lowFrequencyPortfolioHoldOnly

lowFrequencyPortfolioAllDay = portfolio_create(dateStart,dateEnd)
portfolio_settings(lowFrequencyPortfolioAllDay,holdingPeriodsOnly=FALSE)
portfolio_addPosition(lowFrequencyPortfolioAllDay,symbol,quantity=as.numeric(lowFrequencyStrategy),time=printTime)
lowFrequencyPortfolioAllDay

plot1<-util_ggplot(util_plot2d(position_quantity(highFrequencyPortfolioHoldOnly,symbol),title="High Frequency Portfolio Strategy",line_size=0.6))
plot2<-util_ggplot(util_plot2d(position_quantity(lowFrequencyPortfolioHoldOnly,symbol),title="Low Frequency Portfolio Strategy",line_size=0.6))
util_multiplot(plot1,plot2,cols=1)

############################################################
# Part 4 - Trading strategy variance
############################################################
util_plot2d(portfolio_variance(highFrequencyPortfolioHoldOnly),title="Variance, daily",legend="HF HoldOnly")+
util_line2d(portfolio_variance(highFrequencyPortfolioAllDay),legend="HF AllDay")+
util_line2d(portfolio_variance(lowFrequencyPortfolioHoldOnly),legend="LF HoldOnly")+
util_line2d(portfolio_variance(lowFrequencyPortfolioAllDay),legend="LF AllDay")

############################################################
# Part 5 - Trading strategy Value-at-Risk
############################################################

util_plot2d(portfolio_VaR(highFrequencyPortfolioHoldOnly,0.95),title="VaR, daily",legend="HF HoldOnly")+
util_line2d(portfolio_VaR(highFrequencyPortfolioAllDay,0.95),legend="HF AllDay")+
util_line2d(portfolio_VaR(lowFrequencyPortfolioHoldOnly,0.95),legend="LF HoldOnly")+
util_line2d(portfolio_VaR(lowFrequencyPortfolioAllDay,0.95),legend="LF AllDay")

############################################################
# Part 6 - Trading strategy expected return
############################################################

util_plot2d(portfolio_expectedReturn(highFrequencyPortfolioHoldOnly),title="Expected Return, daily",legend="HF HoldOnly")+
util_line2d(portfolio_expectedReturn(highFrequencyPortfolioAllDay),legend="HF AllDay")+
util_line2d(portfolio_expectedReturn(lowFrequencyPortfolioHoldOnly),legend="LF HoldOnly")+
util_line2d(portfolio_expectedReturn(lowFrequencyPortfolioAllDay),legend="LF AllDay")

############################################################
# Part 7 - Trading strategy Sharpe ratio
############################################################

util_plot2d(portfolio_sharpeRatio(highFrequencyPortfolioHoldOnly),title="Sharpe Ratio, daily",legend="HF HoldOnly")+
util_line2d(portfolio_sharpeRatio(highFrequencyPortfolioAllDay),legend="HF AllDay")+
util_line2d(portfolio_sharpeRatio(lowFrequencyPortfolioHoldOnly),legend="LF HoldOnly")+
util_line2d(portfolio_sharpeRatio(lowFrequencyPortfolioAllDay),legend="LF AllDay")
