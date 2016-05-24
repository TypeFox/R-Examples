Fmpd <-
function(dist.mat, my.sample){

Fmpd.sub = function(my.sub.sample){ 
		
## Get the names of species present in the 
## community
		com.names = names(my.sub.sample[my.sub.sample > 0])
		
## Calculate mpd by extracting the lower triangle 
## of a phylogenetic distance matrix comprised of 
## only the species in our community.
		mean(as.dist(dist.mat[com.names, com.names]))

	}


apply(my.sample, MARGIN = 1, Fmpd.sub)

}
