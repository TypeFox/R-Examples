data(lombard)

# univariate change point analysis
bcp.model <- bcp(lombard, burnin=500, mcmc=5000, return.mcmc=TRUE)

# linear model change point analysis 
bcpr.model <- bcp(lombard, cbind(1:100), burnin=500, mcmc=5000, return.mcmc=TRUE)

plot(bcp.model, main="Lombard Milling Data")
plot(bcpr.model, main="Lombard Milling Data (with Regression Model)")
