#  File R/coef.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
coef.stergm <- function(object, ...){list(formation=object$formation.fit$coef,
                                          dissolution=object$dissolution.fit$coef)}
coefficients.stergm <- coef.stergm
