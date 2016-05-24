risk.attribution <-
function(tickers, weights = rep(1,length(tickers)), start, end, data, CompanyList = NULL)
{
  MCTR = mctr(tickers, weights, start, end, data)
  Weight = weights/sum(weights)
  CCTR = Weight*MCTR
  sigma = sum(CCTR)
  CCTR_percent = 100*CCTR/sigma
  Volatility = volatility(tickers, start, end, data)
  output = cbind(Weight, MCTR, CCTR, CCTR_percent, Volatility)
  Portfolio = c(1, NA, sigma, 100, sigma)
  output = data.frame(rbind(output, Portfolio))
  
  if(is.null(CompanyList))
    colnames(output) = c("Weight","MCTR","CCTR","CCTR(%)","Volatility")
  else
  {
    cnames = c(CompanyList[tickers,1],"")
    output = data.frame(cbind(I(cnames),output))
    colnames(output) = c("Company Name","Weight","MCTR","CCTR","CCTR(%)","Volatility")
  }
  
  rownames(output)[nrow(output)] =  "Portfolio"
  return(output)
}