#  File degreenet/R/ayulemle.r
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
pause <- function(){readline(prompt="Pause. Press <Enter> to continue...");invisible()}
# Simulate a Yule distribution over 100
# observations with PDf exponent of 3.5

set.seed(1)
s4 <- simyule(n=100, rho=3.5)
table(s4)
pause()

#
# Calculate the MLE and an asymptotic confidence
# interval for the parameters
#

s4est <- ayulemle(s4)
s4est

pause()
#
# Compute the AICC and BIC for the model
#

llyuleall(v=s4est$theta,x=s4)
