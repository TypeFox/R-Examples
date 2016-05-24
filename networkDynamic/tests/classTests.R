#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2013 the statnet development team
######################################################################

# tests to determine if the networkDynamic class is being set appropriately tested for

library(networkDynamic)
require(testthat)
#check if activate.edge sets class.
net <- network.initialize(5)
net[1,2] <-1
nD <-activate.edges(net,onset=0,terminus=1)
if (!is.networkDynamic(nD)){
  stop("activate edges did not set networkDynamic class")
}
nd <- NULL

net <- network.initialize(5)
net[1,2] <-1
nD <-activate.vertices(net,onset=0,terminus=1)
if (!is.networkDynamic(nD)){
  stop("activate vertices did not set networkDynamic class")
}

expect_true(is.networkDynamic(as.networkDynamic(network.initialize(0))))

# ----- as.network.networkDynamic ------
# test removing class
net <- network.initialize(5)
activate.vertices(net,onset=1,terminus=2)
# is networkDynamic class removed
expect_false(is.networkDynamic(as.network(net)))
# is network class retained?
expect_is(net,'network')

# ----- as.networkDynamic.network ----
expect_is(as.networkDynamic(network.initialize(2)),'networkDynamic')