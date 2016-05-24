#  File degreenet/R/pln.r
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
library(dnet)

# Simulate a Poisson Lognormal distribution over 100
# observations with lognormal mean of -1 and lognormal variance of 1
# This leads to a mean of 1

set.seed(1)
s4 <- simpln(n=1000, v=c(-1,1))
table(s4)

#
# Calculate the MLE and an asymptotic confidence
# interval for the parameters
#

s4est <- aplnmle(s4)
s4est

# Calculate the MLE and an asymptotic confidence
# interval for rho under the Yule model
#

s4yuleest <- ayulemle(s4)
s4yuleest

# Calculate the MLE and an asymptotic confidence
# interval for rho under the Waring model
#

s4warest <- awarmle(s4)
s4warest

#
# Compare the AICC and BIC for the three models
#

llplnall(v=s4est$theta,x=s4)
llyuleall(v=s4yuleest$theta,x=s4)
llwarall(v=s4warest$theta,x=s4)

#
# Show a plot, even though it can be misleading
#
x <- 1:100
plot(x=x,y=dpln(x=x,v=s4est$theta,cutoff=1),type='l',log='xy',
     xlab="degree", ylab="PMF", main="Fit of Poisson-LogNormal")
points(x=x,y=epmf)
