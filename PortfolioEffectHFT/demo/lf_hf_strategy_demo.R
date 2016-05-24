############################################################
# Part 1 - Define trading signals and construct portfolios
############################################################

require(PortfolioEffectHFT)

symbol = "GOOG"
dateStart = "2014-10-13 09:30:00"
dateEnd = "2014-10-14 16:00:00"

# Create function of moving average
MA=function(x,order){
  result<-x
  x1<-c(0,x)
  result[(order):NROW(x)]<-(cumsum(x1)[-(1:(order))]-cumsum(x1)[-((NROW(x1)-order+1):NROW(x1))])/order
  result[1:(order-1)]<-cumsum(x[1:(order-1)])/(1:(order-1))
  return(result-0.0000000001)
}

highFrequencyPortfolio=portfolio_create(fromTime=dateStart,toTime=dateEnd)
lowFrequencyportfolio=portfolio_create(fromTime=dateStart,toTime=dateEnd)

portfolio_addPosition(highFrequencyPortfolio,symbol,1)
price=position_price(highFrequencyPortfolio,symbol)
printTime=price[,1]

highFrequencyStrategy=array(0,dim=NROW(price))
highFrequencyStrategy[price[,"value"]>MA(price[,"value"],150)]<-100
lowFrequencyStrategy=array(0,dim=NROW(price))
lowFrequencyStrategy[price[,"value"]>MA(price[,"value"],800)]<-100

# Add position GOOG to portfolios
portfolio_addPosition(portfolio=highFrequencyPortfolio,symbol=symbol,quantity=highFrequencyStrategy,time=printTime)
portfolio_addPosition(lowFrequencyportfolio,symbol=symbol,quantity=lowFrequencyStrategy,time=printTime)

# Display general information about the portfolio at the end of a dataset
print(highFrequencyPortfolio)
print(lowFrequencyportfolio)
plot(lowFrequencyportfolio)

############################################################
# Part 2 - Holding intervals visualization
############################################################
 
plot1<-util_ggplot(util_plot2d(position_quantity(highFrequencyPortfolio,symbol),title="High Frequency Portfolio Strategy",line_size=0.6))
plot2<-util_ggplot(util_plot2d(position_quantity(lowFrequencyportfolio,symbol),title="Low Frequency Portfolio Strategy",line_size=0.6))
util_multiplot(plot1,plot2,cols=1)

############################################################
# Part 3 - Trading strategy variance
############################################################

util_plot2d(portfolio_variance(highFrequencyPortfolio),title="Variance, daily",legend="HF Portfolio")+
util_line2d(portfolio_variance(lowFrequencyportfolio),legend="LF Portfolio")

############################################################
# Part 4 - Trading strategy Value-at-Risk (daily, 95% c.i.)
############################################################

util_plot2d(portfolio_VaR(highFrequencyPortfolio,0.95),title="Value at Risk in %, daily (95% c.i.)",legend="HF Portfolio")+
util_line2d(portfolio_VaR(lowFrequencyportfolio,0.95),legend="LF Portfolio")

############################################################
# Part 5 - Trading strategy Sharpe ratio (daily)
############################################################

util_plot2d(portfolio_sharpeRatio(highFrequencyPortfolio),title="Sharpe Ratio, daily",legend="HF Portfolio")+
util_line2d(portfolio_sharpeRatio(lowFrequencyportfolio),legend="LF Portfolio")
