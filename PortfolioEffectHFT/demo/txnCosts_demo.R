############################################################
# Part 1 - Define trading signals and construct portfolios
############################################################

require(PortfolioEffectHFT)

# Create function of moving average
MA<-function(x,order){
  result<-x
  x1<-c(0,x)
  result[(order):NROW(x)]<-(cumsum(x1)[-(1:(order))]-cumsum(x1)[-((NROW(x1)-order+1):NROW(x1))])/order
  result[1:(order-1)]<-cumsum(x[1:(order-1)])/(1:(order-1))
  return(result-0.0000000001)
}

portfolio<-portfolio_create("2014-10-10 09:30:00","2014-10-10 16:00:00")

portfolio_addPosition(portfolio,"GOOG",1)
goog<-position_price(portfolio,"GOOG")
printTime<-goog[,1]

# Create two strategy by moving average differents length
highFrequencyStrategy<-array(0,dim=NROW(goog))
highFrequencyStrategy[goog[,"value"]>MA(goog[,"value"],150)]<-100
lowFrequencyStrategy<-array(0,dim=NROW(goog))
lowFrequencyStrategy[goog[,"value"]>MA(goog[,"value"],800)]<-100


numberOfTransaction<-NULL
numberOfTransaction<-rbind(numberOfTransaction,data.frame(number=NROW(diff(highFrequencyStrategy)[diff(highFrequencyStrategy)!=0]),Legends="High Frequency Strategy"))
numberOfTransaction<-rbind(numberOfTransaction,data.frame(number=NROW(diff(lowFrequencyStrategy)[diff(lowFrequencyStrategy)!=0]),Legends="Low Frequency Strategy"))

xlabel<-""
ylabel<-"Number of transaction"
ggplot(numberOfTransaction, aes(x =Legends,y = number))+geom_bar(stat = "identity",fill="#014d64",width=.5)+xlab(xlabel)+ylab(ylabel)+  util_plotTheme()+util_fillScheme()+ggtitle("Number of Transaction")

############################################################
# Part 2 - Holding Intervals Visualization
############################################################

highFrequencyPortfolioWithTransactionCosts<-portfolio_create("SPY","2014-10-10 09:30:00","2014-10-10 16:00:00")
portfolio_settings(highFrequencyPortfolioWithTransactionCosts,txnCostPerShare=0.02)
portfolio_addPosition(highFrequencyPortfolioWithTransactionCosts,"GOOG",quantity=as.numeric(highFrequencyStrategy),time=printTime)
highFrequencyPortfolioWithTransactionCosts

highFrequencyPortfolioWithoutTransactionCosts<-portfolio_create("SPY","2014-10-10 09:30:00","2014-10-10 16:00:00")
portfolio_addPosition(highFrequencyPortfolioWithoutTransactionCosts,"GOOG",quantity=as.numeric(highFrequencyStrategy),time=printTime)
highFrequencyPortfolioWithoutTransactionCosts

lowFrequencyPortfolioWithTransactionCosts<-portfolio_create("SPY","2014-10-10 09:30:00","2014-10-10 16:00:00")
portfolio_settings(lowFrequencyPortfolioWithTransactionCosts,txnCostPerShare=0.02)
portfolio_addPosition(lowFrequencyPortfolioWithTransactionCosts,"GOOG",quantity=as.numeric(lowFrequencyStrategy),time=printTime)
lowFrequencyPortfolioWithTransactionCosts

lowFrequencyPortfolioWithoutTransactionCosts<-portfolio_create("SPY","2014-10-10 09:30:00","2014-10-10 16:00:00")
portfolio_addPosition(lowFrequencyPortfolioWithoutTransactionCosts,"GOOG",quantity=as.numeric(lowFrequencyStrategy),time=printTime)
lowFrequencyPortfolioWithoutTransactionCosts

# Plot portfolio strategy buy and sell changing over time 
plot1<-util_ggplot(util_plot2d(position_quantity(highFrequencyPortfolioWithTransactionCosts,"GOOG"),title="High Frequency Portfolio Strategy",line_size=0.6))
plot2<-util_ggplot(util_plot2d(position_quantity(lowFrequencyPortfolioWithTransactionCosts,"GOOG"),title="Low Frequency Portfolio Strategy",line_size=0.6))
 util_multiplot(plot1,plot2,cols=1)

############################################################
# Part 3 - Trading strategy variance
############################################################

# portfolio and position variance over time
PV1=portfolio_variance(highFrequencyPortfolioWithoutTransactionCosts)
PV2=portfolio_variance(lowFrequencyPortfolioWithoutTransactionCosts)
util_plot2d(portfolio_variance(highFrequencyPortfolioWithTransactionCosts),title="Variance, daily",legend="HF With Transaction Costs")+
  util_line2d(PV1+cbind(array(0,dim=NROW(PV1)),array(1/700000,dim=NROW(PV1))),legend="HF Without Transaction Costs")+
  util_line2d(portfolio_variance(lowFrequencyPortfolioWithTransactionCosts),legend="LF With Transaction Costs")+
  util_line2d(PV2+cbind(array(0,dim=NROW(PV2)),array(1/700000,dim=NROW(PV2))),legend="LF Without Transaction Costs")

############################################################
# Part 4 - Trading strategy expected return
############################################################

# portfolio and position return over time
util_plot2d(portfolio_expectedReturn(highFrequencyPortfolioWithTransactionCosts),title="Expected Return, daily",legend="HF With Transaction Costs")+
  util_line2d(portfolio_expectedReturn(highFrequencyPortfolioWithoutTransactionCosts),legend="HF Without Transaction Costs")+
  util_line2d(portfolio_expectedReturn(lowFrequencyPortfolioWithTransactionCosts),legend="LF With Transaction Costs")+
  util_line2d(portfolio_expectedReturn(lowFrequencyPortfolioWithoutTransactionCosts),legend="LF Without Transaction Costs")

############################################################
# Part 5 - Trading strategy Transactional Costs
############################################################

# portfolio and position Sharpe Ratio over time
txnCosts=portfolio_txnCosts(lowFrequencyPortfolioWithoutTransactionCosts);
util_plot2d(portfolio_txnCosts(highFrequencyPortfolioWithTransactionCosts),title="Transactional Costs",legend="HF With Transaction Costs")+
  util_line2d(portfolio_txnCosts(highFrequencyPortfolioWithoutTransactionCosts),legend="HF Without Transaction Costs")+
  util_line2d(portfolio_txnCosts(lowFrequencyPortfolioWithTransactionCosts),legend="LF With Transaction Costs")+
  util_line2d(txnCosts+cbind(array(0,dim=NROW(txnCosts)),array(5,dim=NROW(txnCosts))),legend="LF Without Transaction Costs")