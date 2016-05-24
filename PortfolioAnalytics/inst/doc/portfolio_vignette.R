### R code from vignette source 'portfolio_vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: portfolio_vignette.Rnw:56-57
###################################################
library(PortfolioAnalytics)


###################################################
### code chunk number 2: portfolio_vignette.Rnw:62-71
###################################################
data(edhec)

# Use the first 4 columns in edhec for a returns object
returns <- edhec[, 1:4]
colnames(returns) <- c("CA", "CTAG", "DS", "EM")
print(head(returns, 5))

# Get a character vector of the fund names
fund.names <- colnames(returns)


###################################################
### code chunk number 3: portfolio_vignette.Rnw:79-83
###################################################
# Specify a portfolio object by passing a character vector for the 
# assets argument.
pspec <- portfolio.spec(assets=fund.names)
print.default(pspec)


###################################################
### code chunk number 4: portfolio_vignette.Rnw:92-97
###################################################
# Add the full investment constraint that specifies the weights must sum to 1.
pspec <- add.constraint(portfolio=pspec, 
                        type="weight_sum", 
                        min_sum=1, 
                        max_sum=1)


###################################################
### code chunk number 5: portfolio_vignette.Rnw:106-116
###################################################
# The full investment constraint can also be specified with type="full_investment"
# pspec <- add.constraint(portfolio=pspec, type="full_investment")

# Another common constraint is that portfolio weights sum to 0.
# This can be specified any of the following ways
# pspec <- add.constraint(portfolio=pspec, type="weight_sum", 
#                         min_sum=0,
#                         max_sum=0)
# pspec <- add.constraint(portfolio=pspec, type="dollar_neutral")
# pspec <- add.constraint(portfolio=pspec, type="active")


###################################################
### code chunk number 6: portfolio_vignette.Rnw:121-137
###################################################
# Add box constraints
pspec <- add.constraint(portfolio=pspec,
                        type="box",
                        min=0.05,
                        max=0.4)

# min and max can also be specified per asset
# pspec <- add.constraint(portfolio=pspec,
#                         type="box",
#                         min=c(0.05, 0, 0.08, 0.1),
#                         max=c(0.4, 0.3, 0.7, 0.55))

# A special case of box constraints is long only where min=0 and max=1
# The default action is long only if min and max are not specified
# pspec <- add.constraint(portfolio=pspec, type="box")
# pspec <- add.constraint(portfolio=pspec, type="long_only")


###################################################
### code chunk number 7: portfolio_vignette.Rnw:143-149
###################################################
# Add group constraints
pspec <- add.constraint(portfolio=pspec, type="group",
                        groups=list(groupA=c(1, 2, 3),
                                    grouB=4),
                        group_min=c(0.1, 0.15), 
                        group_max=c(0.85, 0.55))


###################################################
### code chunk number 8: portfolio_vignette.Rnw:155-160
###################################################
# Add position limit constraint such that we have a maximum number of three assets with non-zero weights.
pspec <- add.constraint(portfolio=pspec, type="position_limit", max_pos=3)

# Can also specify maximum number of long positions and short positions
# pspec <- add.constraint(portfolio=pspec, type="position_limit", max_pos_long=3, max_pos_short=3)


###################################################
### code chunk number 9: portfolio_vignette.Rnw:165-166
###################################################
pspec <- add.constraint(portfolio=pspec, type="diversification", div_target=0.7)


###################################################
### code chunk number 10: portfolio_vignette.Rnw:171-172
###################################################
pspec <- add.constraint(portfolio=pspec, type="turnover", turnover_target=0.2)


###################################################
### code chunk number 11: portfolio_vignette.Rnw:177-178
###################################################
pspec <- add.constraint(portfolio=pspec, type="return", return_target=0.007)


###################################################
### code chunk number 12: portfolio_vignette.Rnw:183-186
###################################################
pspec <- add.constraint(portfolio=pspec, type="factor_exposure",
                        B=c(-0.08, 0.37, 0.79, 1.43),
                        lower=0.6, upper=0.9)


###################################################
### code chunk number 13: portfolio_vignette.Rnw:191-192
###################################################
pspec <- add.constraint(portfolio=pspec, type="transaction_cost", ptc=0.01)


###################################################
### code chunk number 14: portfolio_vignette.Rnw:196-197
###################################################
print(pspec)


###################################################
### code chunk number 15: portfolio_vignette.Rnw:201-202
###################################################
summary(pspec)


###################################################
### code chunk number 16: portfolio_vignette.Rnw:210-243
###################################################
# full investment constraint
weight_constr <- weight_sum_constraint(min_sum=1, max_sum=1)

# box constraint
box_constr <- box_constraint(assets=pspec$assets, min=0, max=1)

# group constraint
group_constr <- group_constraint(assets=pspec$assets, 
                                 groups=list(c(1, 2, 3),
                                             4),
                                 group_min=c(0.1, 0.15), 
                                 group_max=c(0.85, 0.55),
                                 group_labels=c("GroupA", "GroupB"))

# position limit constraint
poslimit_constr <- position_limit_constraint(assets=pspec$assets, max_pos=3)

# diversification constraint
div_constr <- diversification_constraint(div_target=0.7)

# turnover constraint
to_constr <- turnover_constraint(turnover_target=0.2)

# target return constraint
ret_constr <- return_constraint(return_target=0.007)

# factor exposure constraint
exp_constr <- factor_exposure_constraint(assets=pspec$assets,
                                         B=c(-0.08, 0.37, 0.79, 1.43),
                                         lower=0.6, upper=0.9)

# transaction cost constraint
ptc_constr <- transaction_cost_constraint(assets=pspec$assets, ptc=0.01)


###################################################
### code chunk number 17: portfolio_vignette.Rnw:252-256
###################################################
pspec <- add.objective(portfolio=pspec,
                       type='risk',
                       name='ETL',
                       arguments=list(p=0.95))


###################################################
### code chunk number 18: portfolio_vignette.Rnw:261-264
###################################################
pspec <- add.objective(portfolio=pspec,
                       type='return',
                       name='mean')


###################################################
### code chunk number 19: portfolio_vignette.Rnw:269-275
###################################################
pspec <- add.objective(portfolio=pspec, type="risk_budget", name="ETL", 
                       arguments=list(p=0.95), max_prisk=0.3)

# for an equal risk contribution portfolio, set min_concentration=TRUE
# pspec <- add.objective(portfolio=pspec, type="risk_budget", name="ETL", 
#                        arguments=list(p=0.95), min_concentration=TRUE)


###################################################
### code chunk number 20: portfolio_vignette.Rnw:291-293
###################################################
pspec <- add.objective(portfolio=pspec, type="weight_concentration", 
                       name="HHI", conc_aversion=0.1)


###################################################
### code chunk number 21: portfolio_vignette.Rnw:297-302
###################################################
pspec <- add.objective(portfolio=pspec, type="weight_concentration", 
                       name="HHI",
                       conc_aversion=c(0.03, 0.06),
                       conc_groups=list(c(1, 2),
                                        c(3, 4)))


###################################################
### code chunk number 22: portfolio_vignette.Rnw:306-307
###################################################
print(pspec)


###################################################
### code chunk number 23: portfolio_vignette.Rnw:311-312
###################################################
summary(pspec)


###################################################
### code chunk number 24: portfolio_vignette.Rnw:330-362
###################################################
R <- edhec[, 1:4]

# set up simple portfolio with leverage and box constraints 
pspec <- portfolio.spec(assets=colnames(R))
pspec <- add.constraint(portfolio=pspec, type="leverage", 
                        min_sum=0.99, max_sum=1.01)
pspec <- add.constraint(portfolio=pspec, type="box", min=0, max=1)

# generate random portfolios using the 3 methods
rp1 <- random_portfolios(portfolio=pspec, permutations=5000, 
                         rp_method='sample')
rp2 <- random_portfolios(portfolio=pspec, permutations=5000, 
                         rp_method='simplex') 
rp3 <- random_portfolios(portfolio=pspec, permutations=5000, 
                         rp_method='grid')

# show feasible portfolios in mean-StdDev space
tmp1.mean <- apply(rp1, 1, function(x) mean(R %*% x))
tmp1.StdDev <- apply(rp1, 1, function(x) StdDev(R=R, weights=x))
tmp2.mean <- apply(rp2, 1, function(x) mean(R %*% x))
tmp2.StdDev <- apply(rp2, 1, function(x) StdDev(R=R, weights=x))
tmp3.mean <- apply(rp3, 1, function(x) mean(R %*% x))
tmp3.StdDev <- apply(rp3, 1, function(x) StdDev(R=R, weights=x))

# plot feasible portfolios 
plot(x=tmp1.StdDev, y=tmp1.mean, col="gray", main="Random Portfolio Methods",
     ylab="mean", xlab="StdDev")
points(x=tmp2.StdDev, y=tmp2.mean, col="red", pch=2)
points(x=tmp3.StdDev, y=tmp3.mean, col="lightgreen", pch=5)
legend("bottomright", legend=c("sample", "simplex", "grid"), 
       col=c("gray", "red", "lightgreen"),
       pch=c(1, 2, 5), bty="n")


###################################################
### code chunk number 25: portfolio_vignette.Rnw:368-378
###################################################
fev <- 0:5
par(mfrow=c(2, 3))
for(i in 1:length(fev)){
  rp <- rp_simplex(portfolio=pspec, permutations=2000, fev=fev[i])
  tmp.mean <- apply(rp, 1, function(x) mean(R %*% x))
  tmp.StdDev <- apply(rp, 1, function(x) StdDev(R=R, weights=x))
  plot(x=tmp.StdDev, y=tmp.mean, main=paste("FEV =", fev[i]),
       ylab="mean", xlab="StdDev", col=rgb(0, 0, 100, 50, maxColorValue=255))
}
par(mfrow=c(1,1))


###################################################
### code chunk number 26: portfolio_vignette.Rnw:384-400
###################################################
par(mfrow=c(1, 2))
# simplex
rp_simplex <- random_portfolios(portfolio=pspec, permutations=2000, 
                                rp_method='simplex')
tmp.mean <- apply(rp_simplex, 1, function(x) mean(R %*% x))
tmp.StdDev <- apply(rp_simplex, 1, function(x) StdDev(R=R, weights=x))
plot(x=tmp.StdDev, y=tmp.mean, main="rp_method=simplex fev=0:5",
     ylab="mean", xlab="StdDev", col=rgb(0, 0, 100, 50, maxColorValue=255))
#sample
rp_sample <- random_portfolios(portfolio=pspec, permutations=2000, 
                               rp_method='sample')
tmp.mean <- apply(rp_sample, 1, function(x) mean(R %*% x))
tmp.StdDev <- apply(rp_sample, 1, function(x) StdDev(R=R, weights=x))
plot(x=tmp.StdDev, y=tmp.mean, main="rp_method=sample",
     ylab="mean", xlab="StdDev", col=rgb(0, 0, 100, 50, maxColorValue=255))
par(mfrow=c(1,1))


###################################################
### code chunk number 27: portfolio_vignette.Rnw:426-441
###################################################
library(DEoptim)
library(ROI)
require(ROI.plugin.glpk)
require(ROI.plugin.quadprog)

data(edhec)
R <- edhec[, 1:6]
colnames(R) <- c("CA", "CTAG", "DS", "EM", "EQMN", "ED")
funds <- colnames(R)

# Create an initial portfolio object with leverage and box constraints
init <- portfolio.spec(assets=funds)
init <- add.constraint(portfolio=init, type="leverage", 
                       min_sum=0.99, max_sum=1.01)
init <- add.constraint(portfolio=init, type="box", min=0.05, max=0.65)


###################################################
### code chunk number 28: portfolio_vignette.Rnw:446-447
###################################################
maxret <- add.objective(portfolio=init, type="return", name="mean")


###################################################
### code chunk number 29: portfolio_vignette.Rnw:451-456
###################################################
opt_maxret <- optimize.portfolio(R=R, portfolio=maxret, 
                                 optimize_method="ROI", 
                                 trace=TRUE)

print(opt_maxret)


###################################################
### code chunk number 30: portfolio_vignette.Rnw:460-463
###################################################
plot(opt_maxret, risk.col="StdDev", return.col="mean", 
     main="Maximum Return Optimization", chart.assets=TRUE,
     xlim=c(0, 0.05), ylim=c(0,0.0085))


###################################################
### code chunk number 31: portfolio_vignette.Rnw:468-469
###################################################
minvar <- add.objective(portfolio=init, type="risk", name="var")


###################################################
### code chunk number 32: portfolio_vignette.Rnw:473-476
###################################################
opt_minvar <- optimize.portfolio(R=R, portfolio=minvar, 
                                 optimize_method="ROI", trace=TRUE)
print(opt_minvar)


###################################################
### code chunk number 33: portfolio_vignette.Rnw:480-483
###################################################
plot(opt_minvar, risk.col="StdDev", return.col="mean", 
     main="Minimum Variance Optimization", chart.assets=TRUE,
     xlim=c(0, 0.05), ylim=c(0,0.0085))


###################################################
### code chunk number 34: portfolio_vignette.Rnw:488-490
###################################################
qu <- add.objective(portfolio=init, type="return", name="mean")
qu <- add.objective(portfolio=qu, type="risk", name="var", risk_aversion=0.25)


###################################################
### code chunk number 35: portfolio_vignette.Rnw:494-498
###################################################
opt_qu <- optimize.portfolio(R=R, portfolio=qu, 
                             optimize_method="ROI", 
                             trace=TRUE)
print(opt_qu)


###################################################
### code chunk number 36: portfolio_vignette.Rnw:501-504
###################################################
plot(opt_qu, risk.col="StdDev", return.col="mean", 
     main="Quadratic Utility Optimization", chart.assets=TRUE,
     xlim=c(0, 0.05), ylim=c(0, 0.0085))


###################################################
### code chunk number 37: portfolio_vignette.Rnw:509-510
###################################################
etl <- add.objective(portfolio=init, type="risk", name="ETL")


###################################################
### code chunk number 38: portfolio_vignette.Rnw:514-518
###################################################
opt_etl <- optimize.portfolio(R=R, portfolio=etl, 
                              optimize_method="ROI", 
                              trace=TRUE)
print(opt_etl)


###################################################
### code chunk number 39: portfolio_vignette.Rnw:521-524
###################################################
plot(opt_etl, risk.col="ES", return.col="mean", 
     main="ETL Optimization", chart.assets=TRUE,
     xlim=c(0, 0.14), ylim=c(0,0.0085))


###################################################
### code chunk number 40: portfolio_vignette.Rnw:529-532
###################################################
meanETL <- add.objective(portfolio=init, type="return", name="mean")
meanETL <- add.objective(portfolio=meanETL, type="risk", name="ETL",
                         arguments=list(p=0.95))


###################################################
### code chunk number 41: portfolio_vignette.Rnw:536-540
###################################################
opt_meanETL <- optimize.portfolio(R=R, portfolio=meanETL, 
                                  optimize_method="random",
                                  trace=TRUE, search_size=2000)
print(opt_meanETL)


###################################################
### code chunk number 42: portfolio_vignette.Rnw:544-547
###################################################
stats_meanETL <- extractStats(opt_meanETL)
dim(stats_meanETL)
head(stats_meanETL)


###################################################
### code chunk number 43: portfolio_vignette.Rnw:551-553
###################################################
plot(opt_meanETL, risk.col="ETL", return.col="mean", 
     main="mean-ETL Optimization", neighbors=25)


###################################################
### code chunk number 44: portfolio_vignette.Rnw:557-560
###################################################
pct_contrib <- ES(R=R, p=0.95, portfolio_method="component", 
                  weights=extractWeights(opt_meanETL))
barplot(pct_contrib$pct_contrib_MES, cex.names=0.8, las=3, col="lightblue")


###################################################
### code chunk number 45: portfolio_vignette.Rnw:567-576
###################################################
# change the box constraints to long only
init$constraints[[2]]$min <- rep(0, 6)
init$constraints[[2]]$max <- rep(1, 6)

rb_meanETL <- add.objective(portfolio=init, type="return", name="mean")
rb_meanETL <- add.objective(portfolio=rb_meanETL, type="risk", name="ETL",
                            arguments=list(p=0.95))
rb_meanETL <- add.objective(portfolio=rb_meanETL, type="risk_budget", 
                            name="ETL", max_prisk=0.4, arguments=list(p=0.95))


###################################################
### code chunk number 46: portfolio_vignette.Rnw:580-585
###################################################
opt_rb_meanETL <- optimize.portfolio(R=R, portfolio=rb_meanETL, 
                                     optimize_method="DEoptim", 
                                     search_size=2000, 
                                     trace=TRUE, traceDE=5)
print(opt_rb_meanETL)


###################################################
### code chunk number 47: portfolio_vignette.Rnw:588-591
###################################################
plot(opt_rb_meanETL, risk.col="ETL", return.col="mean", 
     main="Risk Budget mean-ETL Optimization",
     xlim=c(0,0.12), ylim=c(0.005,0.009))


###################################################
### code chunk number 48: portfolio_vignette.Rnw:595-597
###################################################
plot.new()
chart.RiskBudget(opt_rb_meanETL, risk.type="percentage", neighbors=25)


###################################################
### code chunk number 49: portfolio_vignette.Rnw:603-609
###################################################
eq_meanETL <- add.objective(portfolio=init, type="return", name="mean")
eq_meanETL <- add.objective(portfolio=eq_meanETL, type="risk", name="ETL",
                            arguments=list(p=0.95))
eq_meanETL <- add.objective(portfolio=eq_meanETL, type="risk_budget", 
                            name="ETL", min_concentration=TRUE, 
                            arguments=list(p=0.95))


###################################################
### code chunk number 50: portfolio_vignette.Rnw:613-618
###################################################
opt_eq_meanETL <- optimize.portfolio(R=R, portfolio=eq_meanETL, 
                                     optimize_method="DEoptim", 
                                     search_size=2000, 
                                     trace=TRUE, traceDE=5)
print(opt_eq_meanETL)


###################################################
### code chunk number 51: portfolio_vignette.Rnw:622-626
###################################################
plot.new()
plot(opt_eq_meanETL, risk.col="ETL", return.col="mean", 
     main="Risk Budget mean-ETL Optimization",
     xlim=c(0,0.12), ylim=c(0.005,0.009))


###################################################
### code chunk number 52: portfolio_vignette.Rnw:630-632
###################################################
plot.new()
chart.RiskBudget(opt_eq_meanETL, risk.type="percentage", neighbors=25)


###################################################
### code chunk number 53: portfolio_vignette.Rnw:644-651
###################################################
opt_combine <- combine.optimizations(list(meanETL=opt_meanETL,
                                          rbmeanETL=opt_rb_meanETL,
                                          eqmeanETL=opt_eq_meanETL))

# View the weights and objective measures of each optimization
extractWeights(opt_combine)
obj_combine <- extractObjectiveMeasures(opt_combine)


###################################################
### code chunk number 54: portfolio_vignette.Rnw:654-655
###################################################
chart.Weights(opt_combine, plot.type="bar", legend.loc="topleft", ylim=c(0, 1))


###################################################
### code chunk number 55: portfolio_vignette.Rnw:659-663
###################################################
plot.new()
chart.RiskReward(opt_combine, risk.col="ETL", return.col="mean", 
                 main="ETL Optimization Comparison", xlim=c(0.018, 0.024),
                 ylim=c(0.005, 0.008))


###################################################
### code chunk number 56: portfolio_vignette.Rnw:667-670
###################################################
STARR <- obj_combine[, "mean"] / obj_combine[, "ETL"]
barplot(STARR, col="blue", cex.names=0.8, cex.axis=0.8,
        las=3, main="STARR", ylim=c(0,1))


###################################################
### code chunk number 57: portfolio_vignette.Rnw:673-676
###################################################
plot.new()
chart.RiskBudget(opt_combine, match.col="ETL", risk.type="percent", 
                 ylim=c(0,1), legend.loc="topright")


