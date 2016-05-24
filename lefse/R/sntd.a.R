sntd.a <-
function(dist.mat, my.sample){

sntd.a.sub = function(x){ 
		
## Get the names of species present in the 
## community.
	com.names = names(x[x > 0])
		
## Make the community phylogenetic distance 
## matrix using the names of the species present ## in the community.
	my.com.dist = dist.mat[com.names, com.names]
		
## Place NA values in the diagonals
	diag(my.com.dist) = NA
		
## Calculate a mean of the minimum values in each 
## row of the community phylogenetic distance 
## matrix weighed by the abundances of the 
## species present in the community.
	wt.sd(apply(my.com.dist, 1, min, na.rm = T), x[x > 0])
	
}

apply(my.sample, MARGIN = 1, sntd.a.sub)

}
