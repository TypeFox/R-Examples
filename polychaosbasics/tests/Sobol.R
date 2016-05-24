###########################################################
# Computation of Sensitivity Indexes (SI)
# for the so-called Sobol function
# via Polynomial Chaos Expansion (PCE)
###########################################################


# Load of necessary functions
library("polychaosbasics")

# nvx<-4 # Number of factors of interest
nlhs<-200 # 20000
degree<-5
#############################################
# Build Legendre polynomial
#############################################
set.seed(42)
pce <- analyticsPolyLeg(nlhs, degree, 'sobol')
# -----------------------------------------
# Default display
# -----------------------------------------
print(pce)

# -----------------------------------------
# Polynomial structure matrix
# -----------------------------------------
print(pce@design@.Data)

#############################################
# PCESI calculation
#############################################
retour <- PCESI(pce)

print(retour, all=TRUE)

# -----------------------------------------
# Individual Monomial SI
# -----------------------------------------
print(retour@IMSI)

# -----------------------------------------
# Regression coefficients
# -----------------------------------------
print(retour@coef)

#############################################
# Plot Y against Y.hat
# and regression line
#############################################
y.hat <-retour@y.hat
y.obs <- pce[, "Y"]
reg <- lm(y.hat ~ y.obs)
jpeg(file="Sobol.jpg")
 plot(y.hat, y.obs,
      xlab="metamodel output", ylab="computer model output",
      main="Sobol test", sub="Scatter plot and regression line")
lines(reg$fitted.values, y.obs)
 dev.off()

