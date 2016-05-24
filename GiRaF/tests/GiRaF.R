# test GiRaF
library(GiRaF)

# Dimension of the lattice
height <- 8
width <- 10

# Interaction parameter
Beta <- 0.6 # Isotropic configuration
# Beta <- c(0.6, 0.6) # Anisotropic configuration when nei = 4
# Beta <- c(0.6, 0.6, 0.6, 0.6) # Anisotropic configuration when nei = 8

# Number of colors
K <- 2 
# Number of neighbors 
G <- 4

# Optional potential on sites
potential <- runif(K,-1,1)
# Optional borders. 
Top <- Bottom <- sample(0:(K-1), width, replace = TRUE)
Left <- Right <- sample(0:(K-1), height, replace = TRUE)
Corner <- sample(0:(K-1), 4, replace = TRUE)

# Partition function for the default setting
NC.mrf(h = height, w = width, param = Beta)

# When specifying the number of colors and neighbors
NC.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta)

# When specifying an optional potential on sites
NC.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta, 
       pot = potential)

# When specifying potential borders. The users will omit to mention all
# the non-existing borders
NC.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta, 
       top = Top, left = Left, bottom = Bottom, right = Right, corner = Corner)

# Exact sampling for the default setting
img <- exact.mrf(h = height, w = width, param = Beta, view = TRUE)

# When specifying the number of colors and neighbors
img <- exact.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta, 
          view = TRUE)

# When specifying an optional potential on sites
img <- exact.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta, 
          pot = potential, view = TRUE)

# When specifying potential borders. The users will omit to mention all
# the non-existing borders
img <- exact.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta, 
          top = Top, left = Left, bottom = Bottom, right = Right, corner = Corner, view = TRUE)

# Algorithm settings
n <- 200
method <- "Gibbs"

# Sampling method for the default setting
img <- sampler.mrf(iter = n, sampler = method, h = height, w = width, 
            param = Beta, view = TRUE)

# Sampling using an existing configuration as starting point
img <- sampler.mrf(iter = n, sampler = method, h = height, w = width, 
            ncolors = K, nei = G, param = Beta, 
            initialise = FALSE, view = TRUE)

# Specifying optional arguments. The users may omit to mention all
# the non-existing borders
img <- sampler.mrf(iter = n, sampler = method, h = height, w = width, 
            ncolors = K, nei = G, param = Beta,
            pot = potential, top = Top, left = Left, bottom = Bottom, 
            right = Right, corner = Corner, view = TRUE)

# Gibbs sampler with sequential updates of the sites. 
img <- sampler.mrf(iter = n, sampler = "Gibbs", h = height, w = width, 
            ncolors = K, nei = G, param = Beta,
            random = FALSE, view = TRUE)
