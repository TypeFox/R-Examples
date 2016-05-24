#  File tests/tergm_parallel.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
library(tergm)

opttest({
  data(florentine)
  net <- flobusiness
  set.seed(1)
  
  mod1 <- stergm(flobusiness, formation= ~edges + degree(3), 
                 dissolution= ~offset(edges),
                 offset.coef.diss=log(9), 
                 targets="formation",
                 estimate="EGMME",
                 control=control.stergm(parallel=2, parallel.type="PSOCK")
  )

}, testname='tergm_parallel')

opttest({
  data(florentine)
  net <- flobusiness
  set.seed(1)
  
  mod1 <- stergm(flobusiness, formation= ~edges + degree(3), 
                 dissolution= ~offset(edges),
                 offset.coef.diss=log(9), 
                 targets="formation",
                 estimate="EGMME",
                 control=control.stergm(parallel=2, parallel.type="MPI")
  )
  
}, testname='tergm_parallel_MPI', testvar="ENABLE_MPI_TESTS")