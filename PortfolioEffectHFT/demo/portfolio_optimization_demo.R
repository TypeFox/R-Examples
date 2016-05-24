############################################################
# Part 1 - Construct buy-and-hold portfolio
############################################################

require(PortfolioEffectHFT)

portfolio<-portfolio_create("SPY", "2014-11-19 09:30:00", "2014-11-19 16:00:00")
portfolio_settings(portfolio,portfolioMetricsMode="price",resultsSamplingInterval='1m')
portfolio_addPosition(portfolio,"AAPL",1000)
portfolio_addPosition(portfolio,"GOOG",1000)
portfolio_addPosition(portfolio,"SPY",1000)


plot(portfolio)

############################################################
# Part 2 - Run portfolio optimization with one constraint
############################################################

optimizer=optimization_goal(portfolio,goal="Variance",direction="minimize")
optimizer=optimization_constraint_portfolioValue(optimizer,10^9)
optimizer=optimization_constraint_expectedReturn(optimizer,constraintType=">=",constraintValue=0)
optimPortfolioOneConstraints=optimization_run(optimizer)

plot(optimPortfolioOneConstraints)

############################################################
# Part 3 - Run portfolio optimization with two constraints
############################################################

optimizer=optimization_constraint_sumOfAbsWeights(optimizer,
                                        constraintType=">=",
                                        constraintValue=0.5,c("AAPL","GOOG"))

optimPortfolioTwoConstraints=optimization_run(optimizer)
plot(optimPortfolioTwoConstraints)

############################################################
# Part 4 - Compute optimal portfolio expected return
############################################################

util_plot2d(portfolio_expectedReturn(optimPortfolioOneConstraints),"Portfolio Expected Return",legend="Optimal Portfolio, with Constraints:\nReturn>=0")+
  util_line2d(portfolio_expectedReturn(optimPortfolioTwoConstraints),legend="Optimal Portfolio, with Constraints:\nReturn>=0, Sum of Abs Weights AAPL and GOOG >=0.5")

############################################################
# Part 5 - Compute optimal portfolio variance
############################################################

util_plot2d(portfolio_variance(optimPortfolioOneConstraints),"Portfolio Variance",legend="Optimal Portfolio, with Constraints:\nReturn>=0")+
  util_line2d(portfolio_variance(optimPortfolioTwoConstraints),legend="Optimal Portfolio, with Constraints:\nReturn>=0, \nSum of Abs Weights AAPL and GOOG >=0.5")

############################################################
# Part 6 - Compute optimal portfolio sum of absolute weights
############################################################

sumOfAbsWeightsOptimPortfolio=abs(position_weight(optimPortfolioOneConstraints,"AAPL")[,2])+abs(position_weight(optimPortfolioOneConstraints,"GOOG")[,2])
sumOfAbsWeightsOptimPortfolioTwoConstraints=abs(position_weight(optimPortfolioTwoConstraints,"AAPL")[,2])+abs(position_weight(optimPortfolioTwoConstraints,"GOOG")[,2])
timeUTC <- position_weight(optimPortfolioOneConstraints,"AAPL")[,1]

util_plot2d(cbind(timeUTC,sumOfAbsWeightsOptimPortfolio),"Portfolio Sum Of Abs Weigth AAPL and GOOG",legend="Optimal Portfolio, with Constraints:\nReturn>=0")+
  util_line2d(cbind(timeUTC,sumOfAbsWeightsOptimPortfolioTwoConstraints),legend="Optimal Portfolio, with Constraints:\nReturn>=0, \nSum of Abs Weights AAPL and GOOG >=0.5")
