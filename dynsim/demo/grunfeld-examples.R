# ---------------------------------------------------------------------------- #
# dynsim Grunfeld examples
# Christopher Gandrud
# 21 July 2015
# ---------------------------------------------------------------------------- #

# ----------------------------- Set Up -------------- ------------------------ #
# Load required packages
library(dynsim)
library(DataCombine)

# Load the Grunfeld data set
data(grunfeld, package = "dynsim")

# Create lagged dependent variable
grunfeld <- slide(grunfeld, Var = "invest", GroupVar = "company", 
                  TimeVar = "year", NewVar = "InvestLag")

# ----------------------------- Estimate and Simulate ------------------------ #
# Estimate linear regression model
M1 <- lm(invest ~ InvestLag + mvalue + kstock, data = grunfeld)

# Create simulation scenarios
attach(grunfeld)
Scen1 <- data.frame(InvestLag = mean(InvestLag, na.rm = TRUE),
                    mvalue = quantile(mvalue, 0.95),
                    kstock = quantile(kstock, 0.95))
Scen2 <- data.frame(InvestLag = mean(InvestLag, na.rm = TRUE),
                    mvalue = mean(mvalue),
                    kstock = mean(kstock))
Scen3 <- data.frame(InvestLag = mean(InvestLag, na.rm = TRUE),
                    mvalue = quantile(mvalue, 0.05),
                    kstock = quantile(kstock, 0.05))
detach(grunfeld)

# Combine into a single list
ScenComb <- list(Scen1, Scen2, Scen3)

# Simulate without shocks
Sim1 <- dynsim(obj = M1, ldv = "InvestLag", scen = ScenComb, n = 20)

# Simulate with shocks

# Keep only the mvalue for the first company for the first 15 years
grunfeldsub <- subset(grunfeld, company == 1)
grunfeldshock <- grunfeldsub[1:15, "mvalue"]

# Create data frame for the shock argument
grunfeldshock <- data.frame(times = 1:15, mvalue = grunfeldshock)

Sim2 <- dynsim(obj = M1, ldv = "InvestLag", scen = ScenComb, n = 15,
               shocks = grunfeldshock)

# Model with interaction
M2 <- lm(invest ~ InvestLag + mvalue*kstock, data = grunfeld)

# Simulate
Sim3 <- dynsim(obj = M2, ldv = "InvestLag", scen = ScenComb, n = 15,
               shocks = grunfeldshock)

# ---------------------------- Plotting -------------------------------------- #
# Basic
dynsimGG(Sim1)

# Specify legend labels
# Create legend labels vector
Labels <- c("95th Percentile", "Mean", "5th Percentile")
dynsimGG(Sim1, leg.name = "Scenarios", leg.labels = Labels, color = "OrRd",
         ylab = "Predicted Real Gross Investment\n")

# Specify colour scale and full labels
dynsimGG(Sim2, leg.name = "Scenarios", leg.labels = Labels, color = "OrRd",
         ylab = "Predicted Real Gross Investment\n", shockplot.var = "mvalue",
         shockplot.ylab = "Firm Value")