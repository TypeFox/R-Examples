# generate artificial species pool / community / trait data. The main purpose is to
# make a simplified artificial data set, to get a feeling for how STEPCAM works

generate.Artificial.Data <- function(numSpecies, numTraits, numCommunities,
occurence_distribution, average_richness, SD_richness, random.Mechanism){
	
  # create species x trait matrix based on settings
	pool_richness <- numSpecies
	n_traits <- numTraits;
	n_communities <- numCommunities
	mechanism_random <- random.Mechanism;	
	
	traits_pool <- matrix(round(rnorm(n = (pool_richness * n_traits), mean = 0, sd = 1), 3),
  ncol = n_traits, nrow = pool_richness)
  
	column1 <- c()
	for(i in 1:pool_richness){
	  column1[i] <- paste("species", i, sep = "")
	}
	row1 <- rep(NA, (n_traits + 1))
	for(i in 1:n_traits + 1){
	  row1[i] <- paste("trait", (i - 1), sep="")
	}
	row1[1] <- "species"
	columns <- matrix(cbind(column1, traits_pool), ncol = (n_traits + 1))
	traits_pool_df <- as.data.frame(columns)
	names(traits_pool_df) <- row1

	comm_spec_matrix <- matrix(rep(0, (n_communities * pool_richness)),
  ncol = pool_richness, nrow = n_communities)

	occurence <- 10^rnorm(pool_richness, mean = 0, sd = occurence_distribution)

	richness_values <- round((rnorm(n_communities, mean = 0, sd = 1) * (SD_richness * pool_richness) + 
  (average_richness * pool_richness)), 0)
	richness_values[which(richness_values > pool_richness)] <- pool_richness
	richness_values[which(richness_values < 0)] <- 0
	for(i in 1:n_communities){
	  if(mechanism_random == T){
	    present_species <- sample(c(1:pool_richness), richness_values[i],
      prob = occurence, replace = F)
	  }else{
	    opt_traits <- rnorm(n_traits, mean = 0, sd = 1)
	    opt_traits_matrix <- c()
	    for(j in 1:n_traits){
	      for(k in 1:pool_richness){
		      l <- (j - 1) * pool_richness + k
		      opt_traits_matrix[l] <- opt_traits[j]
	      }
	    }
	    opt_traits_matrix <- matrix(opt_traits_matrix, ncol = n_traits)
	    traits_distances <- abs(traits_pool - opt_traits_matrix)
	    trait_distances_onedimen <- c()
	    for(m in 1:pool_richness){
	      trait_distances_onedimen[m] <- sqrt(sum(traits_distances[m, ] * traits_distances[m, ]))
	    }
	    inv_distances <- -(trait_distances_onedimen - max(trait_distances_onedimen)) + 10^-99
	    present_species <- sample(c(1:pool_richness), richness_values[i],
      prob = inv_distances, replace = F)
    }
	  comm_spec_matrix[i, c(present_species)] <- 1
	}
	comm_spec_matrix
	if(length(which(colSums(comm_spec_matrix) == 0)) != 0) traits_pool_df <-
  traits_pool_df[-c(which(colSums(comm_spec_matrix) == 0)), ]
	row1 <- c(column1[-c(which(colSums(comm_spec_matrix) == 0))])
	if(length(which(colSums(comm_spec_matrix) == 0)) != 0) comm_spec_matrix <-
  comm_spec_matrix[, -c(which(colSums(comm_spec_matrix) == 0))]
		
	comm_spec_matrix <- as.data.frame(cbind(comm_spec_matrix, rep("", n_communities)))
	names(comm_spec_matrix) <- row1

	plots <- as.data.frame(matrix(as.numeric(as.matrix(comm_spec_matrix[, -ncol(comm_spec_matrix)])),
  ncol = (ncol(comm_spec_matrix) - 1)))

	names(plots) <- traits_pool_df[, 1]

	traits_pool_df[ ,c(2:ncol(traits_pool_df))] <- as.matrix(traits_pool_df[ ,c(2:ncol(traits_pool_df))])
	for(i in 2:ncol(traits_pool_df)){                                                                             
    traits_pool_df[, i] <- as.numeric(traits_pool_df[, i])
	}
	output <- list( traits = traits_pool_df, abundances = plots);
	
	return(output)
}

