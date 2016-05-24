############################################################
# Part 1 - Construct portfolio for optimization
############################################################

require(PortfolioEffectHFT)

portfolio=portfolio_create(fromTime="2014-04-13 9:30:01",
			   toTime="2014-04-16 16:00:00")
portfolio_settings(portfolio,portfolioMetricsMode="price",windowLength = '360m')
portfolio_addPosition(portfolio,'GOOG',1)
portfolio_addPosition(portfolio,'AAPL',1)
portfolio_addPosition(portfolio,'C',1)


plot(portfolio)

############################################################
# Part 2 - Compute theoretical efficient frontier
############################################################

portfolio_settings(portfolio,
		   portfolioMetricsMode="price",
		   windowLength = '360m',
		   resultsSamplingInterval='last')

resultLintner=data.frame(Variance=0,ExpectedReturn=0)

for(x in seq(0,0.015,0.005)){
  optimizer=optimization_goal(portfolio,"Variance","minimize")    
  optimizer=optimization_constraint_portfolioValue(optimizer,10^9)
  optimizer=optimization_constraint_expectedReturn(optimizer,"=",x)
  optimPortfolio=optimization_run(optimizer)
  resultLintner=rbind(resultLintner,c(portfolio_variance(optimPortfolio)[2],portfolio_expectedReturn(optimPortfolio)[2]))
}

resultLintner=resultLintner[-1,]

resultLintner=data.frame(Variance=spline(resultLintner$Variance, n=100)$y,
                         ExpectedReturn=spline(resultLintner$ExpectedReturn, n=100)$y)


ggplot()+geom_path(data=resultLintner, aes(x=Variance,y=ExpectedReturn),size=1.2)+util_plotTheme()+ggtitle("Efficient Frontier")+ylab("Expected Return")

##############################################################
# Part 3 - Compute efficient frontiers of realistic portfolios
##############################################################

resultLintner6000Portfolio=data.frame(Variance=0,ExpectedReturn=0)

for(x in seq(0.004,0.016,0.004)){
  optimizer=optimization_goal(portfolio,"Variance","minimize")
  optimizer=optimization_constraint_portfolioValue(optimizer,6000)
  optimizer=optimization_constraint_expectedReturn(optimizer,"=",x)
  optimPortfolio=optimization_run(optimizer)
  resultLintner6000Portfolio=rbind(resultLintner6000Portfolio,c(portfolio_variance(optimPortfolio)[2],portfolio_expectedReturn(optimPortfolio)[2]))
}

resultLintner6000Portfolio=resultLintner6000Portfolio[-1,]

resultLintner20000Portfolio=data.frame(Variance=0,ExpectedReturn=0)

for(x in seq(0,0.016,0.004)){
  optimizer=optimization_goal(portfolio,"Variance","minimize")
  optimizer=optimization_constraint_portfolioValue(optimizer,20000)
  optimizer=optimization_constraint_expectedReturn(optimizer,"=",x)
  optimPortfolio=optimization_run(optimizer)
  resultLintner20000Portfolio=rbind(resultLintner20000Portfolio,c(portfolio_variance(optimPortfolio)[2],portfolio_expectedReturn(optimPortfolio)[2]))
}

resultLintner20000Portfolio=resultLintner20000Portfolio[-1,]

resultLintner6000Portfolio$legend="$6000 Portfolio"
resultLintner20000Portfolio$legend="$20000 Portfolio"
resultLintner$legend="Theoretical Portfolio"
result=rbind(resultLintner6000Portfolio,resultLintner20000Portfolio,resultLintner)

ggplot()+geom_path(data=result, aes(x=Variance,y=ExpectedReturn,col=legend),size=1.2)+
util_plotTheme()+ggtitle("Efficient Frontier of Theoretical/$20000/$6000  portfolio")+
ylab("Expected Return")+util_colorScheme()

##############################################################
# Part 4 - Compare Markowitz and Lintner efficient frontiers
##############################################################

portfolio_settings(portfolio,
		   portfolioMetricsMode="price",
		   windowLength = '360m',
		   resultsSamplingInterval='last',
		   shortSalesMode = 'markowitz')

resultMarkowitz=data.frame(Variance=0,ExpectedReturn=0)

for(x in seq(0,0.015,0.005)){
  optimizer=optimization_goal(portfolio,"Variance","minimize")
  optimizer=optimization_constraint_portfolioValue(optimizer,10^9)
  optimizer=optimization_constraint_expectedReturn(optimizer,"=",x)
  optimPortfolio=optimization_run(optimizer)
  resultMarkowitz=rbind(resultMarkowitz,c(portfolio_variance(optimPortfolio)[2],portfolio_expectedReturn(optimPortfolio)[2]))
}

resultMarkowitz=resultMarkowitz[-1,]

resultMarkowitz=data.frame(Variance=spline(resultMarkowitz$Variance, n=100)$y,
                           ExpectedReturn=spline(resultMarkowitz$ExpectedReturn, n=100)$y)

resultMarkowitz$legend="Markowitz"
resultLintner$legend="Lintner"
result=rbind(resultMarkowitz,resultLintner)

ggplot()+geom_path(data=result, aes(x=Variance,y=ExpectedReturn,col=legend),size=1.2)+
util_plotTheme()+ggtitle("Markowitz and Lintner Efficient Frontier")+ylab("Expected Return")+
util_colorScheme()