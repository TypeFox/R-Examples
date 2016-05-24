# TODO: Add comment
# 
# Author: cws
###############################################################################

network_as_matrix <- function(network){
	.Call( "network_as_matrix", network, PACKAGE = "NetSim" )
}

as.matrix.NetSimNetwork <- function(x, ...){
	return(network_as_matrix(x))
}

create_network <- function(matrix, directed = TRUE, reflexive = FALSE){
	.Call( "create_network", matrix, directed, reflexive, PACKAGE = "NetSim" )
}

set_tie <- function(network, i, j, value){
	.Call("set_tie", network, i, j, value, PACKAGE = "NetSim")
}

#RcppExport SEXP add_random_ties_to_network(SEXP network, SEXP probability);
add_random_ties_to_network <- function(network, probability = 0.5){
	.Call("add_random_ties_to_network", network, probability, PACKAGE = "NetSim")
}

# RcppExport SEXP add_ring_lattice_to_network(SEXP network, SEXP nReciprocalTies);
add_ring_lattice_to_network <- function(network, nReciprocalTies = 2){
	if (nReciprocalTies <= 0) stop("At least one reciprocal tie per actor needs to be set.")
	.Call("add_ring_lattice_to_network", network, nReciprocalTies, PACKAGE = "NetSim")
}

print.NetSimNetwork <- function(x, ...){
	print("Pointer to NetSim network object:", quote = FALSE)
	print(network_as_matrix(x))
}

