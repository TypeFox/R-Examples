#####
# File:     Chapter09.R
#
# Purpose:  Reproduce Examples in Chapter 9 of the book:
#
#           Millard, SP. (2013).  EnvStats: An R Package for 
#             Environmental Statistics.  Springer-Verlag.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  2013/03/05
#####


#######################################################################################

library(EnvStats)

############################################################

# 9.3 Monte Carlo Simulation 
#---------------------------

# 9.3.1 Simulating the Distribution of the Sum of Two Normal Random Variables
#----------------------------------------------------------------------------

df.100 <- data.frame(simulateMvMatrix(n = 100, seed = 20))
y.100 <- with(df.100, Var.1 + Var.2)

df.10000 <- data.frame(simulateMvMatrix(n = 10000, seed = 20))
y.10000 <- with(df.10000, Var.1 + Var.2)

# Table 9.1
#----------
summaryStats(y.100)

quantile(y.100, probs = c(0.05, 0.95))

summaryStats(y.10000)

quantile(y.10000, probs = c(0.05, 0.95))


# Figure 9.1
#-----------
windows()
hist(y.100, freq = FALSE, col = "cyan", 
	xlab = expression(paste("Y = ", X[1], " + ", X[2])), 
	ylab = "Relative Frequency", 
	main = "Empirical and Theoretical Distributions\nBased on 100 Monte Carlo Trials")
pdfPlot(param.list = list(mean = 0, sd = sqrt(2)), add = TRUE, pdf.lwd = 3)


# Figure 9.2
#-----------
windows()
hist(y.10000, freq = FALSE, breaks = 75, col = "cyan", 
	xlab = expression(paste("Y = ", X[1], " + ", X[2])), 
	ylab = "Relative Frequency", 
	main = "Empirical and Theoretical Distributions\nBased on 10000 Monte Carlo Trials")
pdfPlot(param.list = list(mean = 0, sd = sqrt(2)), add = TRUE, pdf.lwd = 3)


# Clean Up
#---------
rm(df.100, y.100, df.10000, y.10000)

#################################################################################

# 9.4.2 Generating Random Numbers from an Arbitrary Distribution
#---------------------------------------------------------------

# Figure 9.3
#-----------
windows()
cdfPlot(main = paste("Example of Inverse Transformation Method", 
	"for Standard Normal Distribution", sep = "\n"))
usr <- par("usr")
arrows(usr[1], 0.8, x1 = qnorm(0.8), y1 = 0.8, col = "red", lwd = 3)
arrows(qnorm(0.8), 0.8, x1 = qnorm(0.8), y1 = usr[3], col = "red", lwd = 3)

# Clean Up
#---------
rm(usr)

##############################################################################

# 9.4.3 Latin Hypercube Sampling
#-------------------------------

# Figure 9.4
#-----------
windows()
pdfPlot(curve.fill = FALSE, 
	main = paste("Four Equal-Probability Intervals", 
	"for a Standard Normal Distribution", sep = "\n"))
usr <- par("usr")
probs <- c(0.25, 0.5, 0.75)
x <- qnorm(probs)
segments(x0 = x, y0 = usr[3], x1 = x, y1 = dnorm(x),
         col = "red", lwd = 3)

# Clean Up
#---------
rm(usr, probs, x)


# Figure 9.5
#-----------
windows()
cdfPlot(main = paste("Four Equal-Probability Intervals", 
	"for a Standard Normal Distribution", sep = "\n"))
usr <- par("usr")
probs <- c(0.25, 0.5, 0.75)
x <- qnorm(probs)
segments(x0 = usr[1], y0 = probs, x1 = x, y1 = probs,
         col = "red", lwd = 3)
segments(x0 = x, y0 = probs, x1 = x, y1 = usr[3],
         col = "red", lwd = 3)

# Clean Up
#---------
rm(usr, probs, x)


# Figure 9.6
#-----------
windows()
cdfPlot(main = paste("Example of Latin Hypercube Sampling", 
	"with Four Equal-Probability Intervals", sep = "\n"))
usr <- par("usr")
probs <- c(0.25, 0.5, 0.75)
x <- qnorm(probs)
segments(x0 = usr[1], y0 = probs, x1 = x, y1 = probs,
         col = "red", lwd = 3)
segments(x0 = x, y0 = probs, x1 = x, y1 = usr[3],
         col = "red", lwd = 3)

probs <- c(0.04, 0.35, 0.70, 0.89)
x <- qnorm(probs)
arrows(x0 = usr[1], y0 = probs, x1 = x, y1 = probs,
         col = "blue", lty = 6, lwd = 3, angle = 25)
arrows(x0 = x, y0 = probs, x1 = x, y1 = usr[3],
         col = "blue", lty = 6, lwd = 3, angle = 25)

# Clean Up
#---------
rm(usr, probs, x)

######################################################################

# 9.4.4 Example of Simple Random Sampling versus Latin Hypercube Sampling
#------------------------------------------------------------------------

# Figure 9.7
#-----------
x.srs <- simulateVector(50, seed = 798)
windows()
hist(x.srs, freq = FALSE, breaks = 15, 
	ylim = c(0, 0.7), col = "cyan", 
	xlab = "50 Random Numbers Based on SRS", 
	ylab = "Relative Frequency", 
	main = "Results of Simple Random Sampling\n(n = 50)")
pdfPlot(add = TRUE, pdf.lwd = 3)

# Clean Up
#---------
rm(x.srs)


# Figure 9.8
#-----------

x.lhs <- simulateVector(50, sample.method="LHS", seed = 798)
windows()
hist(x.lhs, freq = FALSE, breaks = 15, 
	ylim = c(0, 0.4), col = "cyan", 
	xlab = "50 Random Numbers Based on LHS", 
	ylab = "Relative Frequency", 
	main = "Results of Latin Hypercube Sampling\n(n = 50)")
pdfPlot(add = TRUE, pdf.lwd = 3)

# Clean Up
#---------
rm(x.lhs)

#########################################################################

# 9.4.6 Generating Correlated Multivariate Random Numbers
#--------------------------------------------------------

Intake.fcn <- function(CF = 25, IR, FI = 1, EF, ED = 30, BW = 70, AT = ED * 365)
{ 
	(CF * IR * FI * EF * ED) / (BW * AT) 
}

cors <- c(0, 0.1, 0.5, 0.9)

Intake.mat <- matrix(as.numeric(NA), nrow = 5000, ncol = 4, 
	dimnames = list(NULL, paste("Cor", cors, sep = ".")))

for(j in 1:4) {
	IR.EF.df <- data.frame(simulateMvMatrix(5000, 
	distributions = c(IR = "lnormAlt", EF = "lnormAlt"), 
	param.list = list(IR = list(mean = 0.16, cv = 0.07/0.16), 
		EF = list(mean = 35.5, cv = 25/35.5)), 
	cor.mat = matrix(c(1, cors[j], cors[j], 1), ncol = 2), 
	sample.method = "LHS", seed = 428))
	Intake.mat[, j] <- with(IR.EF.df, Intake.fcn(IR = IR, EF = EF))
}

Results.mat <- 70 * apply(Intake.mat, 2, function(x) 
	c(Mean = mean(x), quantile(x, probs = c(0.5, 0.95, 0.975))))

round(Results.mat, 2)


# Clean up
#---------
rm(Intake.fcn, cors, Intake.mat, j, Results.mat)


#########################################################################

# 9.6 Risk Assessment
#--------------------

# 9.6.3 Example:  Quantifying Variability and Parameter Uncertainty
#------------------------------------------------------------------

# Figure 9.10
#------------
Risk.fcn <- function(C = 500, IR, CF = 1e-6, EF, ED, BW = 70, AT = 25550, CSF = 0.1) {
	CSF * (C * IR * CF * EF * ED) / (BW * AT)
}

IR.vec <- c(50, 100, 200)
ED.max <- c(26, 33, 40)

Risk.mat <- matrix(as.numeric(NA), nrow = 10000, ncol = 3)

set.seed(398)
for(j in 1:3) {
	EF <- simulateVector(n = 10000, distribution = "tri", 
		param.list = list(min = 200, mode = 250, max = 350), 
		sample.method = "LHS")
	ED <- simulateVector(n = 10000, distribution = "lnormTruncAlt", 
		param.list = list(mean = 9, cv = 10/9, min = 0, max = ED.max[j]),
		sample.method = "LHS")
	Risk.mat[, j] <- Risk.fcn(IR = IR.vec[j], EF = EF, ED = ED)
}

ecdfPlot(log10(Risk.mat[, 1]), xlim = c(-7, -4), 
	ecdf.col = "green", xlab = "Risk", xaxt = "n", 
	main = "Empirical CDFs of Risk for 3 Scenarios")
axis(1, at = -7:-4, labels = paste("1e", -7:-4, sep = ""))
ecdfPlot(log10(Risk.mat[, 2]), ecdf.col = "blue", add = TRUE)
ecdfPlot(log10(Risk.mat[, 3]), ecdf.col = "red", add = TRUE)
usr <- par("usr")
legend(x = usr[1], y = 0.9, legend = paste("Case", 1:3), 
	col = c("green", "blue", "red"), lty = 1, lwd = 3, bty = "n")

abline(h = 0.95, lty = 3)
text(x = -6.5, y = 0.97, "95th %'ile")

bounds <- apply(Risk.mat[, c(1, 3)], 2, quantile, probs = 0.95)
log10.bounds <- log10(bounds)

segments(x0 = log10.bounds, x1 = log10.bounds, y0 = 0, y1 = 0.95, lty = 2)
text(x = log10.bounds, y = usr[3]/2, signif(bounds, 2))
arrows(x0 = log10.bounds[1], x1 = log10.bounds[2], y0 = 0.2, y1 = 0.2, 
	code = 3)
text(x = mean(log10.bounds), y = 0.1, "Range of\nUncertainty")


# Figure 9.11
#-----------
Risk.mat.4 <- matrix(as.numeric(NA), nrow = 2000, ncol = 250)

ED.max <- simulateVector(250, distribution = "unif", 
	param.list = list(min = 26, max = 40), 
	sample.method = "LHS", seed = 322, sort = TRUE)

for(j in 1:250) {
	IR <- simulateVector(2000, distribution = "unif", 
		param.list = list(min = 50, max = 200),
		sample.method = "LHS")
	EF <- simulateVector(n = 2000, distribution = "tri", 
		param.list = list(min = 200, mode = 250, max = 350), 
		sample.method = "LHS")
	ED <- simulateVector(n = 2000, distribution = "lnormTruncAlt", 
		param.list = list(mean = 9, cv = 10/9, min = 0, max = ED.max[j]),
		sample.method = "LHS")

	Risk.mat.4[, j] <- Risk.fcn(IR = IR, EF = EF, ED = ED)
}

ecdfPlot(log10(Risk.mat.4[, 1]), xlim = c(-7, -4), 
	ecdf.col = 1, xlab = "Risk", xaxt = "n", 
	main = "Empirical CDFs of Risk for Case 4")
axis(1, at = -7:-4, labels = paste("1e", -7:-4, sep = ""))
for(j in 2:250) {
	ecdfPlot(log10(Risk.mat.4[, j]), ecdf.col = j, add = TRUE)
}


# Figure 9.12
#------------





# Clean Up
#---------
rm(Risk.fcn, ED.max, Risk.mat, Risk.mat.4, j, IR, EF, ED,  
	bounds, log10.bounds, usr)

