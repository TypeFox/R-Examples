#####
# File:     Chapter06.R
#
# Purpose:  Reproduce Examples in Chapter 6 of the book:
#
#           Millard, SP. (2013).  EnvStats: An R Package for 
#             Environmental Statistics.  Springer-Verlag.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  2013/02/17
#####


#######################################################################################

library(EnvStats)

############################################################

# 6.2 Prediction Intervals 
#-------------------------

# Figure 6.1
#-----------
n <- 10
k <- 1
N <- 100
PI.mat <- matrix(nrow = N, ncol = 3)
dimnames(PI.mat) <- list(NULL, c("LPL", "UPL", "Future.Obs"))
set.seed(429)
for(i in 1:N) {
	x <- rnorm(n = 10, mean = 5, sd = 1)
	PI.mat[i, c("LPL", "UPL")] <- predIntNorm(x, conf.level = 0.8)$interval$limits
	PI.mat[i, "Future.Obs"] <- rnorm(1, mean = 5, sd = 1)
}

windows()
plot(1:N, PI.mat[, "LPL"], ylim = range(pretty(range(PI.mat[, c("LPL", "UPL")]))), 
	type = "n", xlab = "Experiment Number", 
	ylab = "Prediction Interval and Future Observation",
	main = paste("Coverage of 80% Prediction Intervals for k=1 Future Observation", 
	"for N(5,1) Distribuiton with n=10", sep = "\n"))
errorBar(x = 1:N,  y = rowMeans(PI.mat[, c("LPL", "UPL")]), 
	lower = PI.mat[, "LPL"], upper = PI.mat[, "UPL"], incr = FALSE, 
	bar.ends = TRUE, gap = FALSE, add = TRUE, 
	bar.ends.size = 0.2, col = "dark grey")
points(1:N, PI.mat[, "Future.Obs"], pch = 18, cex = 0.8)

sum(PI.mat[, "LPL"] <= PI.mat[, "Future.Obs"] &  PI.mat[, "Future.Obs"] <= PI.mat[, "UPL"])

rm(n, k, N, PI.mat, x, i)


###################################################################

# 6.2.1 Prediction Intervals for a Normal Distribution
#-----------------------------------------------------

#
# Prediction Intervals for Future Observations
#
EPA.09.Ex.18.1.arsenic.df

attach(EPA.09.Ex.18.1.arsenic.df)

predIntNorm(Arsenic.ppb[Sampling.Period == "Background"], k = 4, 
		pi.type = "upper", conf.level = 0.95, method = "exact")

predIntNorm(Arsenic.ppb[Sampling.Period == "Background"], k = 4, pi.type = "upper",
	conf.level = 0.95)$interval$limits["UPL"]

detach("EPA.09.Ex.18.1.arsenic.df")


#
# Prediction Intervals for Future Means
#
EPA.09.Ex.18.2.chrysene.df

attach(EPA.09.Ex.18.2.chrysene.df)

predIntNorm(log(Chrysene.ppb)[Well.type == "Background"], n.mean = 4, k = 1, 
	pi.type = "upper", conf.level = 0.99)

mean(log(Chrysene.ppb)[Well.type == "Compliance"])

detach("EPA.09.Ex.18.2.chrysene.df")

###################################################################

# 6.2.2 Prediction Intervals for a Lognormal Distribution
#--------------------------------------------------------

#
# Prediction Intervals for Future Observations
#
attach(EPA.94b.tccb.df)

predIntLnorm(TcCB[Area=="Reference"], k=77, method="exact", pi.type="upper", conf.level=0.95)

sum(TcCB[Area=="Cleanup"] > 2.68)

detach("EPA.94b.tccb.df")


#
# Prediction Intervals for Future Means
#
attach(EPA.09.Ex.18.2.chrysene.df)

predIntLnorm(Chrysene.ppb[Well.type == "Background"], n.geomean = 4, k = 1, 
	pi.type = "upper", conf.level = 0.99)

geomean(Chrysene.ppb[Well.type == "Compliance"])

detach("EPA.09.Ex.18.2.chrysene.df")


############################################################################

# 6.2.3 Prediction Intervals for a Gamma Distribution
#----------------------------------------------------

#
# Prediction Intervals for Future Observations
#
attach(EPA.94b.tccb.df)

predIntGamma(TcCB[Area=="Reference"], k=77, method="exact", 
	pi.type="upper", conf.level=0.95)

detach("EPA.94b.tccb.df")


#
# Prediction Intervals for Future Means
#
attach(EPA.09.Ex.18.2.chrysene.df)

predInt.list <- predIntGamma(Chrysene.ppb[Well.type == "Background"], 
	n.transmean = 4, k = 1, pi.type = "upper", conf.level = 0.99)

predInt.list

trans.power <- predInt.list$interval$normal.transform.power

trans.power

mean.of.trans <- mean(Chrysene.ppb[Well.type == "Compliance"] ^ trans.power)

mean.of.trans ^ (1 / trans.power)


detach("EPA.09.Ex.18.2.chrysene.df")
rm(predInt.list, trans.power, mean.of.trans)

############################################################################

# 6.2.4 Nonparametric Prediction Intervals
#-----------------------------------------

#
# Prediction Intervals for Future Observations
#

# Table 6.2
n <- c(seq(5, 20, by = 5), seq(25, 100, by = 25))
round(100 * predIntNparConfLevel(n = n, m = 3, pi.type = "upper"))

rm(n)

# Prediction limit for the next m=4 trichloroethylene concentrations
# at the compliance well
#-------------------------------------------------------------------
EPA.09.Ex.18.3.TCE.df

with(EPA.09.Ex.18.3.TCE.df, 
	predIntNpar(TCE.ppb[Well.type=="Background"], m = 4, lb = 0, pi.type="upper"))


# Prediction Intervals for Future Medians
#----------------------------------------
EPA.09.Ex.18.4.xylene.df

with(EPA.09.Ex.18.4.xylene.df, 
	predIntNpar(Xylene.ppb[Well.type=="Background"], k = 2, m = 3, 
		lb = 0, pi.type="upper"))

####################################################################################

# 6.3.1 Simultaneous Prediction Intervals for a Normal Distribution
#---------------------------------------------------------------------

nw <- 20
nc <- 10

# Set r = Number of Evaluations per year = 2
r <- 2


# Set Individual Test Type I Error to 
# 1 - (1 - SWFPR)^(1 / (Number of Constituents * Number of Wells))
#
# which translates to setting the confidence limit to 
# (1 - SWFPR)^(1 / (Number of Constituents * Number of Wells))

conf.level <- (1 - 0.1)^(1 / (nc * nw))
conf.level


#
# Prediction intervals based on each of the rules
#

attach(EPA.09.Ex.18.1.arsenic.df)
As.Bkgrd <- Arsenic.ppb[Sampling.Period == "Background"]


# 1-of-2 Rule
predIntNormSimultaneous(As.Bkgrd, k = 1, m = 2, r = r, 
	rule = "k.of.m", pi.type = "upper", 
	conf.level = conf.level)

# 1-of-3 Rule
predIntNormSimultaneous(As.Bkgrd, k = 1, m = 3, r = r, 
	rule = "k.of.m", pi.type = "upper", 
	conf.level = conf.level)$interval$limits["UPL"]

# Modified California Rule
predIntNormSimultaneous(As.Bkgrd, r = r, rule = "Modified.CA", 
	pi.type = "upper", conf.level = conf.level)$interval$limits["UPL"]

# 1-of-2 Rule based on Means of order 2
predIntNormSimultaneous(As.Bkgrd, n.mean = 2, k = 1, m = 2, r = r, 
	rule = "k.of.m", pi.type = "upper", 
	conf.level = conf.level)$interval$limits["UPL"]


#
# Create data frame showing the upper prediction limits for each of the rules, 
# along with the power of detecting a change in concentration of 
# three standard deviations at any of the 20 compliance wells 
# during the course of a year
#

n <- sum(!is.na(As.Bkgrd))
rule.vec <- c("k.of.m", "k.of.m", "Modified.CA", "k.of.m")
n.mean.vec <- c(1, 1, 1, 2)
m.vec <- c(2, 3, 4, 2)
n.rules <- length(rule.vec)
UPL.vec <- rep(as.numeric(NA), n.rules)

for(i in 1:n.rules) {
	UPL.vec[i] <- predIntNormSimultaneous(As.Bkgrd, n.mean = n.mean.vec[i], 
		k = 1, m = m.vec[i], r = r, rule = rule.vec[i], 
		pi.type = "upper", conf.level = conf.level)$interval$limits["UPL"]
}

Power.vec <- predIntNormSimultaneousTestPower(n = n, k = 1, 
		m = m.vec, n.mean = n.mean.vec, r = 2, rule = rule.vec, 
		delta.over.sigma = 3, pi.type = "upper", conf.level = conf.level)
 
data.frame(Rule = rule.vec, k = rep(1, n.rules), 
	m = m.vec, N.Mean = n.mean.vec, UPL = round(UPL.vec, 1), 
	Power = round(Power.vec, 2), 
	Total.Samples = n.mean.vec * m.vec * r)

rm(nw, nc, r, conf.level, As.Bkgrd, n, rule.vec, n.mean.vec, m.vec, 
	n.rules, UPL.vec, Power.vec, i)

detach("EPA.09.Ex.18.1.arsenic.df")


####################################################################################

# 6.3.2 Simultaneous Prediction Intervals for a Lognormal Distribution
#---------------------------------------------------------------------

# Look at backround data:
EPA.09.Ex.19.1.sulfate.df

# Set confidence level:
nw <- 50
nc <- 10

conf.level <- (1 - 0.1)^(1 / (nc * nw))
conf.level


# Compare power of detecting a change in concentration of three standard deviations 
# (on the log scale) at any of the 50 compliance wells during the course of a year 
# for the 1-of-2 rule versus the 1-of-3 rule
predIntNormSimultaneousTestPower(n = 25, k = 1, m = 2:3, r = 2, 
	rule = "k.of.m", delta.over.sigma = 3, pi.type = "upper", 
	conf.level = conf.level)


# Compute upper prediction limit based on 1-of-3 rule
#----------------------------------------------------
with(EPA.09.Ex.19.1.sulfate.df, 
	predIntLnormSimultaneous(Sulfate.mg.per.l, k = 1, m = 3, r = 2, 
		rule = "k.of.m", pi.type = "upper", conf.level = conf.level))

#####################################################################################

# 6.3.3 Simultaneous Prediction Intervals for a Gamma Distribution
#-----------------------------------------------------------------
 
# Compute upper prediction limit based on 1-of-3 rule
#----------------------------------------------------
with(EPA.09.Ex.19.1.sulfate.df, 
	predIntGammaSimultaneous(x = Sulfate.mg.per.l, k = 1, m = 3, r = 2, 
		rule = "k.of.m", pi.type = "upper", 
		conf.level = conf.level))$interval$limits["UPL"]


# Clean up
rm(nw, nc, conf.level)

#####################################################################################

# 6.3.4 Simultaneous Nonparametric Prediction Intervals
#------------------------------------------------------

# Look at the backround and compliance well data:
#------------------------------------------------
EPA.09.Ex.19.5.mercury.df

# n  = 20 based on 4 Background Wells with 5 observations per well
# nw = 10 Compliance Wells
# nc =  5 Constituents
# ne =  1 Evaluation per year

n  <- 20
nw <- 10
nc <-  5
ne <-  1


# Set confidence level and alpha level:
conf.level <- (1 - 0.1)^(1 / nc)
conf.level

alpha <- 1 - conf.level
alpha


# Set the number of future sampling occasions r to the product of 
# the number of compliance wells and the number of evaluations per year
r <- nw * ne


# Create data objects that store information about the 
# six different sampling plans
rule.vec <- c(rep("k.of.m", 3), "Modified.CA", rep("k.of.m", 2))
k.vec    <- rep(1, 6)
m.vec    <- c(2:4, 4, 1, 2)
n.median.vec <- c(rep(1, 4), rep(3, 2))
n.plans <- length(rule.vec)


# Alpha-levels based on using Maximum background observation
alpha.vec.Max <- 1 - predIntNparSimultaneousConfLevel(
	n = n, n.median = n.median.vec, 
	k = k.vec, m = m.vec, r = r, 
	rule = rule.vec, pi.type="upper")


# Alpha-levels based on using 2nd largest background observation
alpha.vec.2nd <- 1 - predIntNparSimultaneousConfLevel(
	n = n, n.median = n.median.vec, 
	k = k.vec, m = m.vec, r = r, 
	rule = rule.vec, pi.type="upper",
	n.plus.one.minus.upl.rank = 2)

# Alpha-levels based on using 3rd largest background observation
alpha.vec.3rd <- 1 - predIntNparSimultaneousConfLevel(
	n = n, n.median = n.median.vec, 
	k = k.vec, m = m.vec, r = r, 
	rule = rule.vec, pi.type="upper",
	n.plus.one.minus.upl.rank = 3)


# Now create a data frame listing all of the plans, their associated 
# Type I error rate, and their associated upper prediction limit

Bkgd.Hg.Sorted <- with(EPA.09.Ex.19.5.mercury.df, 
	sort(Mercury.ppb[Well.type == "Background"], 
		decreasing = TRUE))

Candidate.Plans.df <- data.frame(Rule = rep(rule.vec, 3), 
	k = rep(k.vec, 3), m = rep(m.vec, 3),
	Median.n = rep(n.median.vec, 3), 
	Order.Statistic = rep(c("Max", "2nd", "3rd"), each = n.plans), 
	Achieved.alpha = round(
		c(alpha.vec.Max, alpha.vec.2nd, alpha.vec.3rd), 4),
	BG.Limit = rep(Bkgd.Hg.Sorted[1:3], each = n.plans))

Candidate.Plans.df


#--------------------------------------------------------------------
# Eliminate plans that do not achieve the per-constituent alpha level
#--------------------------------------------------------------------
index <- Candidate.Plans.df$Achieved.alpha <= alpha
Candidate.Plans.df <- Candidate.Plans.df[index, ]

Candidate.Plans.df


#-------------------------------------------------
# Determine whether Well CW-1 is out of compliance
#-------------------------------------------------

EPA.09.Ex.19.5.mercury.df[Well == "CW-1", c("Well", "Event", "Mercury.ppb.orig")]
#   Well Event Mercury.ppb.orig
#25 CW-1     1             0.22
#26 CW-1     2              0.2
#27 CW-1     3              <.2
#28 CW-1     4             0.25
#29 CW-1     5             0.24
#30 CW-1     6              <.2


#        Rule Median.n k m Order.Statistic Achieved.alpha BG.Limit Result
# ----------- -------- - - --------------- -------------- -------- ------
#      k.of.m        1 1 3             Max         0.0055     0.28 Pass. First value of 0.22 < BG Limit.
#      k.of.m        1 1 4             Max         0.0009     0.28 Pass. First value of 0.22 < BG Limit.
# Modified.CA        1 1 4             Max         0.0140     0.28 Pass. First value of 0.22 < BG Limit.
#      k.of.m        3 1 2             Max         0.0060     0.28 Pass. First value of 0.22 < BG Limit.
#      k.of.m        1 1 4             2nd         0.0046     0.25 Pass. First value of 0.22 < BG Limit.
#      k.of.m        1 1 4             3rd         0.0135     0.24 Pass. First value of 0.22 < BG Limit.




#-------------------------------------------------
# Determine whether Well CW-2 is out of compliance
#-------------------------------------------------

EPA.09.Ex.19.5.mercury.df[Well == "CW-2", c("Well", "Event", "Mercury.ppb.orig")]
#   Well Event Mercury.ppb.orig
#31 CW-2     1             0.36
#32 CW-2     2             0.41
#33 CW-2     3             0.28
#34 CW-2     4             0.45
#35 CW-2     5             0.43
#36 CW-2     6             0.54


# Resampling is required for all plans because initial value of 0.36
# is greater than BG Limit, and median of first 3 values (also 0.36)
# is greater than BG Limit.

#        Rule Median.n k m Order.Statistic Achieved.alpha BG.Limit Result
# ----------- -------- - - --------------- -------------- -------- ------
#      k.of.m        1 1 3             Max         0.0055     0.28 Pass. Third value of 0.28 = BG Limit.
#      k.of.m        1 1 4             Max         0.0009     0.28 Pass. Third value of 0.28 = BG Limit.
# Modified.CA        1 1 4             Max         0.0140     0.28 Fail. Only 1 of observations 2-4 <= 0.28.
#      k.of.m        3 1 2             Max         0.0060     0.28 Fail. Both medians (0.36 and 0.45) > BG Limit.
#      k.of.m        1 1 4             2nd         0.0046     0.25 Fail. Observations 1-4 all > 0.25.
#      k.of.m        1 1 4             3rd         0.0135     0.24 Fail. Observations 1-4 all > 0.25.



#------------------------------------------------------------------------------
# Assume a normal distribution and compute the power of detecting a change in 
# concentration of three standard deviations at any of the 10 compliance wells 
# during the course of a year, as well as the total number of potential samples 
# that may have to be taken.
#------------------------------------------------------------------------------

Power.vec <- predIntNparSimultaneousTestPower(n = n,
	n.median = Candidate.Plans.df[, "Median.n"],
	k = Candidate.Plans.df[, "k"], m = Candidate.Plans.df[, "m"],
	r = r, rule = as.character(Candidate.Plans.df[, "Rule"]), 
	n.plus.one.minus.upl.rank = 
		match(Candidate.Plans.df[, "Order.Statistic"] , 
			c("Max", "2nd", "3rd")), 
	delta.over.sigma = 3, pi.type = "upper", r.shifted = 1, 
	distribution = "norm", method = "approx")


data.frame(Candidate.Plans.df[, 
	c("Rule", "k", "m", "Median.n", "Order.Statistic")], 
	Power = round(Power.vec, 2), 
	Total.Samples = Candidate.Plans.df$Median.n * Candidate.Plans.df$m * ne)


rm(n, nw, nc, ne, conf.level, alpha, r, rule.vec, k.vec, m.vec, 
	n.median.vec, n.plans, alpha.vec.Max, alpha.vec.2nd, alpha.vec.3rd,
	Bkgd.Hg.Sorted, Candidate.Plans.df, index, Power.vec)


############################################################

# 6.4 Tolerance Intervals 
#------------------------


# 6.4.1 Tolerance Intervals for a Normal Distribution
#----------------------------------------------------

# First assume facility is in compliance/assessment monitoring
# so compute a lower 99% confidence limit for the 95th percentile 
# for the distribution at each of the three compliance wells.
# This is equivalent to computing a lower beta-content tolerance 
# limit with coverage 5% and associated confidence level of 99%.  
#----------------------------------------------------------------

Aldicarb <- EPA.09.Ex.21.1.aldicarb.df$Aldicarb.ppb
Well <- EPA.09.Ex.21.1.aldicarb.df$Well

tolIntNorm(Aldicarb[Well == "Well.1"], coverage = 0.05, 
	ti.type = "lower", conf.level = 0.99)

tolIntNorm(Aldicarb[Well == "Well.2"], coverage = 0.05, 
	ti.type = "lower", 
	conf.level = 0.99)$interval$limits["LTL"]

tolIntNorm(Aldicarb[Well == "Well.3"], coverage = 0.05, 
	ti.type = "lower", 
	conf.level = 0.99)$interval$limits["LTL"]

# One command to compute the LTL for all three wells at once:
sapply(split(Aldicarb, Well), function(x) {
	tolIntNorm(x, coverage = 0.05, ti.type = "lower", 
	conf.level = 0.99)$interval$limits["LTL"]})

# Now assume assume the site is in corrective action monitoring, so 
# compute the one-sided 99% upper confidence limit for the 95th 
# percentile and compare that to the MCL.  This is equivalent to 
# computing an upper beta-content tolerance limit with coverage 
# 95% and associated confidence level of 99%.
#-------------------------------------------------------------------

sapply(split(Aldicarb, Well), function(x) {
	tolIntNorm(x, coverage = 0.95, ti.type = "upper", 
	conf.level = 0.99)$interval$limits["UTL"]})

rm(Aldicarb, Well)

########################################################################

# 6.4.2 Tolerance Intervals for a Lognormal Distribution
#-------------------------------------------------------

# "Hot-Measurment Comparison" for TcCB data
#------------------------------------------
with(EPA.94b.tccb.df, 
	tolIntLnorm(TcCB[Area == "Reference"], coverage = 0.95, ti.type = "upper", 
		conf.level = 0.95))



# Compare chrysene data from three compliance wells
# to upper tolerance limit based on data from two 
# background wells
#--------------------------------------------------

EPA.09.Ex.17.3.chrysene.df

with(EPA.09.Ex.17.3.chrysene.df, 
	tolIntLnorm(Chrysene.ppb[Well.type == "Background"], 
		coverage = 0.95, ti.type = "upper", 
		conf.level = 0.95))$interval$limits["UTL"]


########################################################################

# 6.4.3 Tolerance Intervals for a Gamma Distribution
#---------------------------------------------------

# "Hot-Measurment Comparison" for TcCB data
#------------------------------------------
with(EPA.94b.tccb.df, 
	tolIntGamma(TcCB[Area == "Reference"], coverage = 0.95, ti.type = "upper", 
		conf.level = 0.95))$interval$limits["UTL"]

########################################################################

# 6.4.4 Nonparametric Tolerance Intervals
#----------------------------------------

EPA.09.Ex.17.4.copper.df

with(EPA.09.Ex.17.4.copper.df, 
	tolIntNpar(Copper.ppb[Well.type == "Background"], 
		conf.level = 0.95, ti.type = "upper", lb = 0))




