### R code from vignette source 'ROI_vignette.Rnw'

###################################################
### code chunk number 1: ROI_vignette.Rnw:26-32
###################################################
suppressMessages(library(PortfolioAnalytics))
suppressMessages(library(foreach))
suppressMessages(library(iterators))
suppressMessages(library(ROI))
suppressMessages(library(ROI.plugin.quadprog))
suppressMessages(library(ROI.plugin.glpk))


###################################################
### code chunk number 2: ROI_vignette.Rnw:37-46
###################################################
data(edhec)

# Use the first 4 columns in edhec for a returns object
returns <- edhec[, 1:4]
colnames(returns) <- c("CA", "CTAG", "DS", "EM")
print(head(returns, 5))

# Get a character vector of the fund names
funds <- colnames(returns)


###################################################
### code chunk number 3: ROI_vignette.Rnw:64-75
###################################################
# Create portfolio object
portf_maxret <- portfolio.spec(assets=funds)

# Add constraints to the portfolio object
portf_maxret <- add.constraint(portfolio=portf_maxret, type="full_investment")
portf_maxret <- add.constraint(portfolio=portf_maxret, type="box",
                               min=c(0.02, 0.05, 0.03, 0.02),
                               max=c(0.55, 0.6, 0.65, 0.5))

# Add objective to the portfolio object
portf_maxret <- add.objective(portfolio=portf_maxret, type="return", name="mean")


###################################################
### code chunk number 4: ROI_vignette.Rnw:79-81
###################################################
print(portf_maxret)
summary(portf_maxret)


###################################################
### code chunk number 5: ROI_vignette.Rnw:86-89
###################################################
# Run the optimization
opt_maxret <- optimize.portfolio(R=returns, portfolio=portf_maxret, 
                                 optimize_method="ROI", trace=TRUE)


###################################################
### code chunk number 6: ROI_vignette.Rnw:93-94
###################################################
print(opt_maxret)


###################################################
### code chunk number 7: ROI_vignette.Rnw:98-99
###################################################
summary(opt_maxret)


###################################################
### code chunk number 8: ROI_vignette.Rnw:104-105
###################################################
names(opt_maxret)


###################################################
### code chunk number 9: ROI_vignette.Rnw:109-110
###################################################
extractStats(opt_maxret)


###################################################
### code chunk number 10: ROI_vignette.Rnw:114-115
###################################################
extractWeights(opt_maxret)


###################################################
### code chunk number 11: ROI_vignette.Rnw:120-121
###################################################
plot(opt_maxret, chart.assets=TRUE, xlim=c(0.02, 0.18))


###################################################
### code chunk number 12: ROI_vignette.Rnw:127-129
###################################################
chart.RiskReward(opt_maxret,return.col="mean", risk.col="sd", 
                 chart.assets=TRUE, xlim=c(0.01, 0.05), main="Maximum Return")


###################################################
### code chunk number 13: ROI_vignette.Rnw:134-138
###################################################
bt_maxret <- optimize.portfolio.rebalancing(R=returns, portfolio=portf_maxret,
                                            optimize_method="ROI", 
                                            rebalance_on="quarters", 
                                            training_period=36)


###################################################
### code chunk number 14: ROI_vignette.Rnw:156-164
###################################################
# Create portfolio object
portf_minvar <- portfolio.spec(assets=funds)

# Add full investment constraint to the portfolio object
portf_minvar <- add.constraint(portfolio=portf_minvar, type="full_investment")

# Add objective to minimize variance
portf_minvar <- add.objective(portfolio=portf_minvar, type="risk", name="var")


###################################################
### code chunk number 15: ROI_vignette.Rnw:170-174
###################################################
# Run the optimization
opt_gmv <- optimize.portfolio(R=returns, portfolio=portf_minvar, 
                              optimize_method="ROI", trace=TRUE)
print(opt_gmv)


###################################################
### code chunk number 16: ROI_vignette.Rnw:178-182
###################################################
bt_gmv <- optimize.portfolio.rebalancing(R=returns, portfolio=portf_minvar,
                                         optimize_method="ROI", 
                                         rebalance_on="quarters", 
                                         training_period=36)


###################################################
### code chunk number 17: ROI_vignette.Rnw:190-202
###################################################
# Add long only constraints
portf_minvar <- add.constraint(portfolio=portf_minvar, type="box", 
                               min=0, max=1)

# Add group constraints
portf_minvar <- add.constraint(portfolio=portf_minvar, 
                               type="group", 
                               groups=list(groupA=1,
                                           groupB=c(2, 3),
                                           groupC=4), 
                               group_min=c(0, 0.25, 0.10), 
                               group_max=c(0.45, 0.6, 0.5))


###################################################
### code chunk number 18: ROI_vignette.Rnw:206-210
###################################################
# Run the optimization
opt_minvar <- optimize.portfolio(R=returns, portfolio=portf_minvar, 
                                 optimize_method="ROI", trace=TRUE)
print(opt_minvar)


###################################################
### code chunk number 19: ROI_vignette.Rnw:214-218
###################################################
bt_minvar <- optimize.portfolio.rebalancing(R=returns, portfolio=portf_minvar,
                                            optimize_method="ROI", 
                                            rebalance_on="quarters", 
                                            training_period=36)


###################################################
### code chunk number 20: ROI_vignette.Rnw:234-255
###################################################
# Create initial portfolio object
init_portf <- portfolio.spec(assets=funds)

# Create full investment constraint
fi_constr <- weight_sum_constraint(type="full_investment")

# Create long only constraint
lo_constr <- box_constraint(type="long_only", assets=init_portf$assets)

# Combine the constraints in a list
qu_constr <- list(fi_constr, lo_constr)

# Create return objective
ret_obj <- return_objective(name="mean")

# Create variance objective specifying a risk_aversion parameter which controls
# how much the variance is penalized
var_obj <- portfolio_risk_objective(name="var", risk_aversion=0.25)

# Combine the objectives into a list
qu_obj <- list(ret_obj, var_obj)


###################################################
### code chunk number 21: ROI_vignette.Rnw:260-266
###################################################
# Run the optimization
opt_qu <- optimize.portfolio(R=returns, portfolio=init_portf, 
                             constraints=qu_constr, 
                             objectives=qu_obj, 
                             optimize_method="ROI",
                             trace=TRUE)


###################################################
### code chunk number 22: ROI_vignette.Rnw:270-276
###################################################
bt_qu <- optimize.portfolio.rebalancing(R=returns, portfolio=init_portf,
                                        constraints=qu_constr, 
                                        objectives=qu_obj, 
                                        optimize_method="ROI", 
                                        rebalance_on="quarters", 
                                        training_period=36)


