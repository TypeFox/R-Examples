############################################################
# Part 1 - Compute optimal weights using Treynor-Black model
############################################################

require(PortfolioEffectHFT)

# create test portfolio
timeStart="2014-10-02 09:30:00"
timeEnd="2014-10-03 16:00:00"

portfolio=portfolio_create("SPY", timeStart, timeEnd)
portfolio_settings(portfolio,portfolioMetricsMode="price",jumpsModel='all')
portfolio_addPosition(portfolio,"AAPL",100)
portfolio_addPosition(portfolio,"GOOG",100)
portfolio_addPosition(portfolio,"SPY",100)

util_plot2d(position_jensensAlpha(portfolio,"AAPL"), title="Jensen's Alpha",legend="AAPL")+
util_line2d(position_jensensAlpha(portfolio,"GOOG"), legend="GOOG")

# compute optimal weights according to the Treynor-Black model
timeUTC=position_jensensAlpha(portfolio,"AAPL")[,1]
alpha=cbind (position_jensensAlpha(portfolio,"AAPL")[,2], position_jensensAlpha(portfolio,"GOOG")[,2])
variance=cbind(position_variance(portfolio,"AAPL")[,2], position_variance(portfolio,"GOOG")[,2])

treynorBlack=alpha/variance
optimWeigth=treynorBlack/rowSums(abs(treynorBlack))

# plot optimal position weights
util_plot2d(cbind(timeUTC[!(is.na(optimWeigth[,1]))], optimWeigth[!(is.na(optimWeigth[,1])),1]), title="Optimal Weight",legend="AAPL")+
util_line2d(cbind(timeUTC[!(is.na(optimWeigth[,2]))], optimWeigth[!(is.na(optimWeigth[,2])),2]), legend="GOOG")

################################################
# Part 2- Compare two portfolios of equal value
################################################

# compute optimal position quatities for a portfolio of given size
portfolioCash=10000000
optimPosition=portfolioCash*optimWeigth/cbind(position_price(portfolio,"AAPL")[,2],position_price(portfolio,"GOOG")[,2])

portfolioSimple=portfolio_create("SPY", timeStart, timeEnd)
portfolio_settings(portfolioSimple,portfolioMetricsMode="price",jumpsModel='all')
portfolio_addPosition(portfolioSimple, "AAPL", quantity = (portfolioCash/2)%/%(position_price(portfolio,"AAPL")[,2]), time = timeUTC)
portfolio_addPosition(portfolioSimple, "GOOG", quantity = (portfolioCash/2)%/%(position_price(portfolio,"GOOG")[,2]), time = timeUTC)

portfolioOptimal=portfolio_create("SPY", timeStart, timeEnd)
portfolio_settings(portfolioOptimal,portfolioMetricsMode="price",jumpsModel='all')
portfolio_addPosition(portfolioOptimal,"AAPL", quantity = optimPosition[,1], time = timeUTC)
portfolio_addPosition(portfolioOptimal,"GOOG", quantity = optimPosition[,2], time = timeUTC)

meanAAPL=mean(position_price(portfolioOptimal,"AAPL")[,2])
meanGOOG=mean(position_price(portfolioOptimal,"GOOG")[,2])

portfolioSimpleAlpha=portfolio_jensensAlpha(portfolioSimple)
portfolioOptimAlpha=portfolio_jensensAlpha(portfolioOptimal)

util_plot2d(portfolioSimpleAlpha, title="Jensen's Alpha",legend="Simple portfolio")+
util_line2d(portfolioOptimAlpha, legend="Optimal portfolio")+
util_line2d(cbind(timeUTC, mean(portfolioSimpleAlpha[,2])), legend="Avg. of simple")+
util_line2d(cbind(timeUTC, mean(portfolioOptimAlpha[,2])), legend="Avg. of optimal")

############################################################
# Part 3 - Use ARMA-class model to forecast alpha-decay
############################################################
require(forecast)
meanSimple=NULL
forecastErrorsSimple=NULL

#forecast test for simple portfolio
#use ARIMA(1,1,0). 1 second foreacast.
for(x1 in seq(0,2,0.1)){
  x2=1-x1
  position_setQuantity(portfolioSimple, "AAPL", (portfolioCash*x1)%/%meanAAPL)
  position_setQuantity(portfolioSimple, "GOOG", (portfolioCash*x2)%/%meanGOOG)
  
  alpha=portfolio_jensensAlpha(portfolioSimple)
  meanSimple=c(meanSimple,mean(alpha[,2]))
  
  forecastErrors=array(0,dim=100)
  for(i in 1:100){
    if((i/21+x1*1000/21)%%5==0){
      print(paste(i/21+x1*1000/21,"%"))
    }
    fit=arima(alpha[(3200+i*400):(3400+i*400-1),2],c(1,1,0)) #ARIMA estimation 
    forecastErrors[i]=abs(forecast(fit,1)$mean-alpha[3400+i*400,2])/mean(abs(alpha[(3200+i*400):(3400+i*400-1),2])) #Calculation of forecast errors 
  } 
  forecastErrorsSimple=c(forecastErrorsSimple,mean(forecastErrors))
}

#forecast test for optimal portfolio
for(i in 1:100){
  if(i%%5==0){
    print(paste(i,"%"))
  }
  fit=arima(portfolioOptimAlpha[(3200+i*400):(3400+i*400-1),2],c(1,1,0)) #ARIMA estimation 
  forecastErrors[i]=abs(forecast(fit,1)$mean-portfolioOptimAlpha[3400+i*400,2])/abs(mean(portfolioOptimAlpha[(3200+i*400):(3400+i*400-1),2])) #Calculation of forecast errors 
} 
forecastErrorsOptim=mean(forecastErrors) #Calculation of average values 

resultdf=data.frame(err=c(forecastErrorsSimple,forecastErrorsOptim) ,mean=c(meanSimple,mean(portfolioOptimAlpha[,2])),legend=c(array("Portfolio Simple Alpha",dim=21),"Portfolio Optimal Alpha"))

ggplot() + geom_point(data=resultdf, aes(x=err, y=mean,colour=legend),size=5) +
  xlab("Forecast Error(%)") + 
  ylab("Alpha mean") +
  util_plotTheme(base_size = 15)+util_colorScheme()
