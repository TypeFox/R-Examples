############################################################
# Part 1 - Construct buy-and-hold portfolio
############################################################

require(PortfolioEffectHFT)


dateStart = "2014-10-01 09:30:00"
dateEnd = "2014-10-02 16:00:00"

# Create empty portfolio
portfolio<-portfolio_create(dateStart,dateEnd)

# Add position AAPL and GOOG to portfolio
portfolio_addPosition(portfolio,"AAPL",100)
portfolio_addPosition(portfolio,"GOOG",200)

# Display general information about the portfolio at the end of the period
print(portfolio)
plot(portfolio)

############################################################
# Part 2 - Compute portfolio value
############################################################

print(util_plot2d(portfolio_value(portfolio),title="Portfolio value, in USD"))

############################################################
# Part 3 - Compute portfolio & position expected return (daily)
############################################################

print(util_plot2d(position_expectedReturn(portfolio,"AAPL"),title="Expected Return, daily",legend="AAPL")+
  util_line2d(position_expectedReturn(portfolio,"GOOG"),legend="GOOG")+
  util_line2d(portfolio_expectedReturn(portfolio),legend="Portfolio"))

############################################################
# Part 4 - Compute portfolio & position variance
############################################################

# Compute portfolio and position variance (daily)
print(util_plot2d(position_variance(portfolio,"AAPL"),title="Variance, daily",legend="AAPL")+
  util_line2d(position_variance(portfolio,"GOOG"),legend="GOOG")+
  util_line2d(portfolio_variance(portfolio),legend="Portfolio"))

############################################################
# Part 5 - Compute portfolio & position Value-at-Risk (daily, 95% c.i.)
############################################################

print(util_plot2d(position_VaR(portfolio,"AAPL",0.95),title="Value at Risk in %, daily (95% c.i.)",legend="AAPL")+
  util_line2d(position_VaR(portfolio,"GOOG",0.95),legend="GOOG")+
  util_line2d(portfolio_VaR(portfolio,0.95),legend="Portfolio"))

############################################################
# Part 6 - Compute portfolio & position Sharpe Ratio (daily)
############################################################

print(util_plot2d(position_sharpeRatio(portfolio,"AAPL"),title="Sharpe Ratio, daily",legend="AAPL")+
  util_line2d(position_sharpeRatio(portfolio,"GOOG"),legend="GOOG")+
  util_line2d(portfolio_sharpeRatio(portfolio),legend="Portfolio"))
