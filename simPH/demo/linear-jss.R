#############
# Linear effect examples from JSS paper
# Christopher Gandrud
# Updated 5 April 2014
#############

# Load packages
library("survival")
library("simPH")
library("ggplot2")
library("gridExtra")

#### Illustration of linear effects####
# Load hmohiv data from UCLA repository
hmohiv <- read.table(
           "http://www.ats.ucla.edu/stat/r/examples/asa/hmohiv.csv", 
           sep = ",", header = TRUE)

# Center age at its median (35)
hmohiv$AgeMed <- hmohiv$age - 35

M1 <- coxph(Surv(time, censor) ~ AgeMed + drug, 
            method = "breslow", data = hmohiv)

# Simulate relative hazards
Sim1 <- coxsimLinear(M1, b = "AgeMed", Xj = seq(-15, 19, by = 1))

# Plot results with simGG default
Plot1_1 <- simGG(Sim1)

# Plot results
Plot1_2 <-simGG(Sim1, xlab = "\nYears of Age from the\n Sample Median (35)",
                ylab = "Relative Hazard with Comparison\n to a 35 Year Old\n",
                alpha = 0.05, type = 'lines')

# Combine plots
grid.arrange(Plot1_1, Plot1_2, ncol = 2)

# Simulate and plot binary drug relative hazard estiamates
Sim2 <- coxsimLinear(M1, b = "drug", Xj = 0:1)

simGG(Sim2, psize = 3, xlab = "",
      ylab = "Relative Hazard\n",
      type = 'points', method = 'lm') + 
    scale_x_continuous(breaks = c(0, 1), 
                       labels = c('\nNo Drug Use', '\nDrug Use'))
