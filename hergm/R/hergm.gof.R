###########################################################################
# Copyright 2009 Nobody                                                   #
#                                                                         #
# This file is part of hergm.                                             #
#                                                                         # 
#    hergm is free software: you can redistribute it and/or modify        #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         # 
#    hergm is distributed in the hope that it will be useful,             #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
#    GNU General Public License for more details.                         #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with hergm.  If not, see <http://www.gnu.org/licenses/>.       #
#                                                                         # 
###########################################################################

hergm.gof <- function(sample = NULL, 
                      verbose = 1,
                      ...)
# input: postprocessed sample, number of nodes, number of blocks, number covariates, observed network
# output: goodness of fit data
{
  # Extract
  network <- sample$network
  n <- sample$n
  model <- sample$model
  max_number <- sample$max_number
  indicator <- sample$indicator
  d1 <- sample$d1
  d2 <- sample$d2
  ergm_theta <- sample$ergm_theta
  hergm_theta <- sample$hergm_theta
  sample_size <- sample$sample_size
  
  # Initialize
  component.number <- vector(length = sample_size) 
  max.component.size <- vector(length = sample_size)
  distance.label <- matrix(0, nrow = sample_size, ncol = n)
  distance <- matrix(0, nrow = sample_size, ncol = n)
  edges <- vector(length = sample_size)
  degree <- matrix(0, nrow = sample_size, ncol = n)
  stars <- vector(length = sample_size)
  triangles <- vector(length = sample_size)

  # Sample
  for (i in 1:sample_size)
    { 
    eta <- c(ergm_theta[i,], hergm_theta[i,]) 
    if (verbose == 1) 
      {
      cat("Sample", i)
      cat("\n------------------------------------------------------------------")
      cat("\nInput:")
      if (d1 > 0) cat("\n- parameters:", ergm_theta[i,])
      if (d2 > 0) 
        {
        cat("\n- block parameters:", hergm_theta[i,])
        cat("\n- block memberships:", indicator[i,])
        }
      cat("\n")
      }
    sample.network <- hergm(  
      model$formula, 
      max_number = max_number, 
      eta = eta, 
      indicator = indicator[i,], 
      simulate = TRUE, 
      samplesize = 1,
      predictions = TRUE,
      verbose = verbose,
    )
    # Construct simulated network:
    simulated_edgelist <- cbind(sample.network$heads, sample.network$tails) # Edge list of simulated network
    if (is.null(simulated_edgelist) == FALSE) 
      {
      simulated_network <- as.network(simulated_edgelist, directed = is.directed(network), matrix.type = "edgelist") # The simulated network as network object; note: simulated_network$gal$n <- maximum vertex number
      if (simulated_network$gal$n < n) add.vertices(simulated_network, nv = n - simulated_network$gal$n) # If simulated$gal$n < n, add isolates 
      # Functions of data:
      components <- component.dist(simulated_network)
      component.number[i] <- length(components$csize) # Number of components
      max.component.size[i] <- max(components$csize) # Size of largest component
      distances <- geodist(simulated_network)
      distances <- distances$gdist
      frequencies_distances <- table(distances) # First column: frequency of self loops; columns 2:number: frequencies finite and (last column) infinite distances
      number <- length(frequencies_distances) - 1 - (sum(distances == Inf) > 0) # Number of distances minus 0-distance minus Inf-distance
      distance.label[i,1:number] <- rownames(frequencies_distances)[2:(number+1)]
      distance[i,1:number] <- frequencies_distances[2:(number+1)] # Frequencies of finite distances
      if (is.directed(network)) # Directed networks
        { 
        edges[i] <- summary(simulated_network ~ edges)
        degree[i,] <- summary(simulated_network ~ odegree(1:n-1)) # Degree distribution          
        stars[i] <- summary(simulated_network ~ ostar(2))     
        triangles[i] <- summary(simulated_network ~ ttriple)
        } 
      else # Undirected networks
        {
        edges[i] <- summary(simulated_network ~ edges)
        degree[i,] <- summary(simulated_network ~ degree(1:n-1)) # Degree distribution          
        stars[i] <- summary(simulated_network ~ kstar(2))     
        triangles[i] <- summary(simulated_network ~ triangles)
        }
      }
    if (verbose) 
      {
      cat("\nnumber of components:", component.number[i])
      cat("\nsize of largest component:", max.component.size[i])
      cat("\nfrequencies of geodesic distances:", distance[i,1:number])
      if (is.directed(network))
        {
        cat("\nnumber of edges:", edges[i])
        cat("\noutdegree distribution:", degree[i,])
        cat("\nnumber of 2-out-stars:", stars[i])
        cat("\nnumber of transitive triples:", triangles[i])
        }
      else 
        {
        cat("\nnumber of edges:", edges[i])
        cat("\ndegree distribution:", degree[i,])
        cat("\nnumber of 2-stars:", stars[i])
        cat("\nnumber of triangles:", triangles[i])
        } 
      cat("\n\n")
      }
    }

  # Store
  gof <- list()
  gof$component.number <- component.number
  gof$max.component.size <- max.component.size
  gof$distance.label <- distance.label
  gof$distance <- distance
  gof$edges <- edges
  gof$degree <- degree
  gof$stars <- stars
  gof$triangles <- triangles

  gof
}

