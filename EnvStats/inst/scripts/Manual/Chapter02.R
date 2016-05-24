#####
# File:     Chapter02.R
#
# Purpose:  Reproduce Examples in Chapter 2 of the book:
#
#           Millard, SP. (2013).  EnvStats: An R Package for 
#             Environmental Statistics.  Springer-Verlag.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  2013/02/03
#####


#######################################################################################

library(EnvStats)

############################################################

# 2.8.1	Confidence Interval for the Mean of a Normal Distribution
#----------------------------------------------------------------
EPA.09.Ex.16.1.sulfate.df

Sulfate.back <- with(EPA.09.Ex.16.1.sulfate.df,
	Sulfate.ppm[Well.type == "Background"])
 
enorm(Sulfate.back, ci = TRUE)

enorm(Sulfate.back, ci = TRUE, ci.param = "variance")$interval

ciNormHalfWidth(n.or.n1 = 8, sigma.hat = c(15, 30, 60))

ciNormHalfWidth(n.or.n1 = 4, sigma.hat = 30)

ciNormN(half.width = 10, sigma.hat = 30)

# Figure 2.1
#-----------
windows()
plotCiNormDesign(sigma.hat=30, range.x.var=c(4, 80), 
  conf=0.99, xlim=c(0, 80), ylim=c(0, 90), 
  main = paste("Half-Width vs. Sample Size for CI for Mean",
    "with Estimated SD=30 and Various Confidence Levels", sep = "\n"))
plotCiNormDesign(sigma.hat=30, range.x.var=c(4, 80), conf=0.95, plot.lty=2, 
	plot.col = 2, add= TRUE)
plotCiNormDesign(sigma.hat=30, range.x.var=c(4, 80), conf=0.90, plot.lty=3, 
	plot.col = 4, add= TRUE)
legend("topright", paste(c("99%", "95%", "90%"), "Confidence"), col= c(1, 2, 4), 
	lty=1:3, lwd=3, bty = "n")

#--------------------------------------------------------------------------------------

Sulfate.down <- with(EPA.09.Ex.16.1.sulfate.df,
	Sulfate.ppm[Well.type == "Downgradient"])

enorm(Sulfate.down)

t.test(Sulfate.down, Sulfate.back, var.equal = TRUE)$conf.int

summary(lm(Sulfate.ppm ~ Well.type, data = EPA.09.Ex.16.1.sulfate.df))$sigma

ciTableMean(n1 = 8, n2 = 8, diff = c(100, 50, 0), 
  SD = c(15, 25, 50), digits = 0)


rm(Sulfate.back, Sulfate.down)

############################################################

# 2.8.2	Confidence Interval for a Binomal Proportion
#---------------------------------------------------

EPA.92c.benzene1.df

with(EPA.92c.benzene1.df, ebinom(Censored, ci = TRUE))

ciBinomHalfWidth(n.or.n1 = 36, p.hat = c(0.75, 0.85, 0.95))

ciBinomHalfWidth(n.or.n1 = 10, p.hat = 0.9)

ciBinomN(half.width = 0.03, p.hat = 0.9)

# Figure 2.2
#-----------
windows()
plotCiBinomDesign(p.hat=0.9, range.x.var=c(10, 200), conf=0.99, 
  xlim=c(0, 200), ylim=c(0, 0.3), 
  main=paste("Half-Width vs. Sample Size for CI for p", 
  "With Estimated p=0.9 and Various Confidence Levels", sep = "\n"))
plotCiBinomDesign(p.hat=0.9, range.x.var=c(10, 200), conf=0.95, plot.col = 2, plot.lty=2, add= TRUE)
plotCiBinomDesign(p.hat=0.9, range.x.var=c(10, 200), conf=0.90, plot.col = 4, plot.lty=3, add= TRUE)
legend("topright", paste(c("99%", "95%", "90%"), "Confidence"), col = c(1, 2, 4), lty=1:3, lwd=3, bty = "n")
mtext("CI Method = Score normal approximation, with continuity correction", line = 0)

#--------------------------------------------------------------------------------------

ciTableProp(n1 = 36, p1.hat = c(0.2, 0.4, 0.6), 
  n2 = 36, p2.hat.minus.p1.hat = c(0.3, 0.15, 0))

############################################################

# 2.8.3	Nonparametric Confidence Interval for a Percentile
#---------------------------------------------------------

EPA.92c.copper2.df

Cu.Bkgrd <- with(EPA.92c.copper2.df, Copper[Well.type == "Background"])

eqnpar(Cu.Bkgrd, p = 0.95)

ciNparConfLevel(n=24, p=0.95, ci.type="upper")

ciNparConfLevel(n=12, p=0.95, ci.type="upper")

ciNparN(p=0.95, ci.type="upper")

# Figure 2.3
#-----------
windows()
plotCiNparDesign(p=0.95, ci.type="upper", 
  range.x.var=c(2, 100), ylim=c(0, 1))


rm(Cu.Bkgrd)

############################################################

# 2.9.1	Prediction Interval for a Normal Distribution
#----------------------------------------------------

EPA.92c.arsenic3.df

As.Bkgrd <- with(EPA.92c.arsenic3.df, Arsenic[Well.type == "Background"])
predIntNorm(As.Bkgrd, k = 4, method = "exact")


# Figure 2.4
#-----------
windows()
plotPredIntNormDesign(range.x.var=c(4, 50), k=4, sigma.hat=17, 
  xlim=c(0, 50), ylim=c(30, 110), 
  main=paste("Half-Width vs. Sample Size for Prediction Interval", 
  "With Estimate SD=17 and Various Numbers of Future Values", sep = "\n"))
plotPredIntNormDesign(range.x.var=c(4, 50), k=2, sigma.hat=17, 
	plot.col = 2, plot.lty=2, add= TRUE)
plotPredIntNormDesign(range.x.var=c(4, 50), k=1, sigma.hat=17, 
	plot.col = 4, plot.lty=3, add= TRUE)
legend("topright", c("k=4", "k=2", "k=1"), col = c(1, 2, 4), lty=1:3, lwd = 3, bty = "n")

rm(As.Bkgrd)

############################################################

# 2.9.2 Nonparametric Prediction Interval
#----------------------------------------

predIntNparN(m = rep(c(1, 5, 10), 2), 
  conf.level = rep(c(0.9, 0.95), each = 3))

# Figure 2.5
#-----------
windows()
plotPredIntNparDesign(range.x.var = c(2, 100), k = 1, m = 1, 
  xlim = c(0, 100), ylim = c(0, 1),
  main=paste("Confidence Level vs. Sample Size for Nonparametric PI", 
  "With Various Values of m (Two-Sided PI)", sep = "\n"))
plotPredIntNparDesign(range.x.var=c(2, 100), k=5, m=5, plot.col = 2, plot.lty=2, add= TRUE)
plotPredIntNparDesign(range.x.var=c(2, 100), k=10, m=10, plot.col = 4, plot.lty=3, add= TRUE)
legend("bottomright", c("m=  1", "m=  5", "m=10"), col = c(1, 2, 4), lty = 1:3, lwd=3, bty = "n")


############################################################

# 2.10.1 Tolerance Interval for a Normal Distribution
#----------------------------------------------------


As.Bkgrd <- with(EPA.92c.arsenic3.df, Arsenic[Well.type == "Background"])
tolIntNorm(As.Bkgrd,  coverage = 0.95, cov.type = "content", 
    ti.type = "two-sided", conf.level = 0.99)

# Figure 2.6
#-----------
windows()
plotTolIntNormDesign(range.x.var = c(5, 50), sigma.hat = 17, 
  coverage = 0.99, conf = 0.99, 
  xlim = c(0, 50), ylim = c(0, 200), 
  main=paste("Half-Width vs. Sample Size for TI With Estimated", 
  "SD = 17 and Various Values of Coverage (99% Confidence)", 
  sep = "\n"))
plotTolIntNormDesign(range.x.var=c(5, 50), sigma.hat=17, coverage=0.95, 
  conf=0.99, plot.col = 2, plot.lty=2, add= TRUE)
plotTolIntNormDesign(range.x.var=c(5, 50), sigma.hat=17, coverage=0.90, 
  conf=0.99, plot.col = 4, plot.lty=3, add= TRUE)
legend("topright", paste(c("99%", "95%", "90%"), "Coverage"), col = c(1, 2, 4), 
	lty=1:3, lwd=3, bty = "n")

rm(As.Bkgrd)

############################################################

# 2.10.2 Nonparametric Tolerance Interval
#----------------------------------------

tolIntNparN(coverage = rep(c(0.8, 0.9, 0.95), 2), 
  conf.level = rep(c(0.9, 0.95), each = 3))


# Figure 2.7
#-----------
windows()
plotTolIntNparDesign(range.x.var = c(2, 100), coverage = 0.8, 
  xlim = c(0, 100), ylim = c(0, 1), 
  main=paste("Confidence Level vs. Sample Size for Nonparametric TI", 
  "With Various Values of Coverage (Two-Sided TI)", 
  sep = "\n"))
plotTolIntNparDesign(range.x.var = c(2, 100), coverage = 0.90, 
  plot.col = 2, plot.lty = 2, add = TRUE)
plotTolIntNparDesign(range.x.var = c(2, 100), coverage = 0.95, 
  plot.col = 4, plot.lty = 3, add = TRUE)
legend("bottomright", paste(c("80%", "90%", "95%"), "Coverage"), 
  col = c(1, 2, 4), lty = 1:3, lwd = 3, bty = "n")

#########################################################################

# 2.11.1 Testing the Mean of a Normal Distribution
#-------------------------------------------------

summaryStats(VC.ppb ~ Well, data = EPA.09.Ex.22.1.VC.df, 
  subset = Period == "Background", digits = 1)

VC.lm.fit <- lm(VC.ppb ~ Well, data = EPA.09.Ex.22.1.VC.df, 
  subset = Period == "Background")

summary(VC.lm.fit)$sigma

sqrt(enorm(VC.lm.fit$residuals, ci = TRUE, ci.param = "variance")$interval$limits)

tTestAlpha(n.or.n1 = 4, delta.over.sigma = 5 / c(3.2, 6), 
  power = 0.8, sample.type = "one.sample", 
  alternative = "greater")

tTestPower(n.or.n1 = c(4, 8, 12), delta.over.sigma = 5 / 3.2, 
  alpha = 0.01, alternative = "greater")

tTestN(delta.over.sigma = 5 / 6, alpha = 0.01, power = 0.9, 
  alternative = "greater")

rm(VC.lm.fit)


# Figure 2.8
#-----------
windows()
plotTTestDesign(alpha = 0.01, delta.over.sigma = 1, range.x.var = c(2, 35), 
  xlim = c(0, 35), ylim = c(0, 1), alternative="greater", approx = FALSE)


# Figure 2.9
#-----------
windows()
plotTTestDesign(y.var = "delta.over.sigma", alpha = 0.01, power = 0.9, 
  range.x.var = c(2, 15), xlim = c(0, 15), ylim = c(0, 40), 
  alternative = "greater", approx = FALSE)


#########################################################################

# 2.11.2 Testing a Binomial Proportion
#-------------------------------------------------

# Sample size based on normal approximation to the binomial
# without continuity correction
#----------------------------------------------------------
propTestN(p.or.p1 = 0.3, p0.or.p2 = 0.1, alpha = 0.05, 
  power = 0.9, sample.type = "one.sample", 
  alternative = "greater", approx = TRUE, round.up = TRUE)

# Simulation showing true Type I error rate using test based on
# normal approximation WITHOUT continuity correction is about 7%, not 5%
#-----------------------------------------------------------------------
set.seed(274)
N <- 10000
Reject.vec <- logical(N)
for(i in 1:N) {
  Reject.vec[i] <- prop.test(
    x = rbinom(n = 1, size = 30, prob = 0.1), 
    n = 30, p = 0.1, alternative = "greater", correct = FALSE)$p.value < 0.05
}
mean(Reject.vec)

binom.test(x = sum(Reject.vec), n = length(Reject.vec), p = 0.05)$conf.int


# Simulation showing true Type I error rate using test based on
# normal approximation WITH continuity correction is about 2.6%, not 5%
#-----------------------------------------------------------------------
set.seed(538)
N <- 10000
Reject.vec <- logical(N)
for(i in 1:N) {
  Reject.vec[i] <- prop.test(
    x = rbinom(n = 1, size = 30, prob = 0.1), 
    n = 30, p = 0.1, alternative = "greater", correct = TRUE)$p.value < 0.05
}
mean(Reject.vec)

binom.test(x = sum(Reject.vec), n = length(Reject.vec), p = 0.05)$conf.int

rm(N, Reject.vec, i)


# Sample size based on exact binomial test
#-----------------------------------------
propTestN(p.or.p1 = 0.3, p0.or.p2 = 0.1, alpha = 0.05, 
  power = 0.9, sample.type = "one.sample", 
  alternative = "greater", approx = FALSE, round.up = TRUE)

propTestN(p.or.p1 = 0.3, p0.or.p2 = 0.1, alpha = 0.05, 
  tol.alpha = 0.01, power = 0.9, sample.type = "one.sample", 
  alternative = "greater", approx = FALSE, round.up = TRUE)


#########################################################################


# 2.11.3 Testing Multiple Wells for Compliance with Simultaneous Prediction Intervals
#------------------------------------------------------------------------------------
nc <- 20
nw <- 100
conf.level <- (1 - 0.1)^(1 / (nc * nw))
conf.level

n <- 25
r <- 2

# Compute power of detecting an increase of three standard deviations 
# in concentration using the prediction interval based on the “1-of-2” 
# resampling rule.
#---------------------------------------------------------------------
predIntNormSimultaneousTestPower(n = n, k = 1, m = 2, r = r, 
  rule = "k.of.m", delta.over.sigma = 3, pi.type = "upper", 
  conf.level = conf.level)


# Reproduce table shown in Step 2 on page 19-23 of the EPA (2009) 
# guidance document.
#----------------------------------------------------------------
rule.vec <- c(rep("k.of.m", 3), "Modified.CA", rep("k.of.m", 3))
m.vec <- c(2, 3, 4, 4, 1, 2, 1)
n.mean.vec <- c(rep(1, 4), 2, 2, 3)
n.scenarios <- length(rule.vec)

K.vec <- predIntNormSimultaneousK(
		n = n, k = 1, m = m.vec,  
		n.mean = n.mean.vec,
		r = r, rule = rule.vec, 
		pi.type = "upper", conf.level = conf.level)

Power.vec <- predIntNormSimultaneousTestPower(
		n = n, k = 1, m = m.vec, 
		n.mean = n.mean.vec, r = r, 
		rule = rule.vec, delta.over.sigma = 3, 
		pi.type = "upper", conf.level = conf.level)

data.frame(Rule = rule.vec, k = rep(1, n.scenarios), 
	m = m.vec, N.Mean = n.mean.vec, 
	K = round(K.vec, 2), 
	Power = round(Power.vec, 2), 
	Total.Samples = m.vec * n.mean.vec)




# Figure 2.10
#------------
# Power for designs based on SINGLE future observations
#------------------------------------------------------
windows()
par(mar = c(5,4,6,1) + 0.1)
plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 4, r = r, rule="k.of.m", pi.type = "upper", 
	conf.level = conf.level,  
	xlab = "SD Units Above Background",
	main = "")
plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 3, r = r, rule="k.of.m", pi.type = "upper", 
	conf.level = conf.level,  
	add = TRUE, plot.col = 2, plot.lty = 2)
plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 2, r = r, rule="k.of.m", pi.type = "upper", 
	conf.level = conf.level,  
	add = TRUE, plot.col = 3, plot.lty = 3)
plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	r = r, rule="Modified.CA", pi.type = "upper", 
	conf.level = conf.level,  
	add = TRUE, plot.col = 4, plot.lty = 4)
legend(0, 1, c("1-of-4", "Modified CA", "1-of-3", "1-of-2"), 
	col = c(1, 4, 2, 3), lty = c(1, 4, 2, 3), lwd = 2, bty="n")
mtext(paste('Power of 1-of-2, 1-of-3, 1-of-4, and Modified CA Designs to',
	"\nDetect An Increase At", nw, "Wells for SWFPR = 10%"), 
	line = 3, cex = 1.25)
mtext(paste("(Based on", n, "Background Observations and\n", 
	"Adjusted for", nc, "Consituents and", r, "Evaluations per year)"), 
	line = 1)


# Figure 2.11
#------------
# Power for designs based on the MEAN of future observations
#-----------------------------------------------------------
windows()
par(mar = c(5,4,6,1) + 0.1)
plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 2, n.mean = 2, r = r, rule="k.of.m", 
	pi.type = "upper", conf.level = conf.level,  
	xlab = "SD Units Above Background",
	main = "")
plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 1, n.mean = 2, r = r, rule="k.of.m", 
	pi.type = "upper", 
	conf.level = conf.level, 
	add = TRUE, plot.col = 2, plot.lty = 2)
plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 1, n.mean = 3, r = r, rule="k.of.m", 
	pi.type = "upper", conf.level = conf.level, 
	add = TRUE, plot.col = 3, plot.lty = 3)
legend(0, 1, c("1-of-2, Order 2", "1-of-1, Order 3", "1-of-1, Order 2"), 
	col = c(1, 3, 2), lty = c(1, 3, 2), lwd = 2, 
	bty="n")
mtext(paste('Power of 1-of-1 and 1-of-2 Designs Based on Means to Detect',
	"\nAn Increase At", nw, "Wells for SWFPR = 10%"), 
	line = 3, cex = 1.25)
mtext(paste("(Based on", n, "Background Observations and\n", 
	"Adjusted for", nc, "Consituents and", r, "Evaluations per Year)"), 
	line = 1)


rm(nc, nw, conf.level, n, r, rule.vec, m.vec, n.mean.vec, 
  n.scenarios, K.vec, Power.vec)



