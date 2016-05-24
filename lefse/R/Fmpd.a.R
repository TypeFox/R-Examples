Fmpd.a <-
function(dist.mat, my.sample){

Fmpd.a.sub = function(x){ 

## Get the names of the species in the community.
	com.names = names(x[x > 0])
		
## Make a matrix with one row containing 
## abundances and names of all present species
	com = t(as.matrix(x[x > 0]))
		
## Make phylogenetic distance matrix for taxa in ## community.
	com.dist = dist.mat[com.names, com.names]
		
## Calculate the product of the abundances of all 
## species in the community.
	abundance.products = t(as.matrix(com[1, com[1, ] > 0, drop = F]))%*% as.matrix(com[1, com[1, ] > 0, drop = F])
		
## Calculate a mean of the community phylogenetic 
## distance matrix weighted by the products of 
## all pairwise abundances.
	weighted.mean(com.dist, abundance.products)
	
	}


apply(my.sample, MARGIN = 1, Fmpd.a.sub)

}
