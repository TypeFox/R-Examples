## ----echo = FALSE-----------------------------------------
library(GiRaF)
library(knitr)
options(prompt = " ", continue='  ', width=60)

## ----eval = FALSE-----------------------------------------
#  # NC.mrf
#  # Dimension of the lattice
#  height <- 8
#  width <- 10
#  
#  # Interaction parameter
#  Beta <- 0.6 # Isotropic configuration
#  # Beta <- c(0.6, 0.6) # Anisotropic configuration for a first
#                        # order dependency structure (nei = 4).
#  # Beta <- c(0.6, 0.6, 0.6, 0.6) # Anisotropic configuration for a second
#                                  # order dependency structure (nei = 8).
#  
#  # Number of colors. Automatically set to 2 if not specified.
#  K <- 2
#  # Number of neighbors. Automatically set to 4 if not specified.
#  G <- 4
#  
#  # Optional potential on sites. Automatically set to NULL if not specified
#  potential <- runif(K,-1,1)
#  # Optional borders. Automatically set to NULL if not specified
#  Top <- Bottom <- sample(0:(K-1), width, replace = TRUE)
#  Left <- Right <- sample(0:(K-1), height, replace = TRUE)
#  Corner <- sample(0:(K-1), 4, replace = TRUE)
#  
#  # Partition function for the default setting
#  NC.mrf(h = height, w = width, param = Beta)
#  # When specifying the number of colors and neighbors
#  NC.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta)
#  # When specifying an optional potential on sites
#  NC.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta,
#         pot = potential)
#  # When specifying possible borders. The users will omit to mention all
#  # the non-existing borders
#  NC.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta,
#         top = Top, left = Left, bottom = Bottom, right = Right,
#         corner = Corner)

## ----eval = FALSE-----------------------------------------
#  # exact.mrf
#  # Dimension of the lattice
#  height <- 8
#  width <- 10
#  
#  # Interaction parameter
#  Beta <- 0.6 # Isotropic configuration
#  # Beta <- c(0.6, 0.6) # Anisotropic configuration for a first
#                        # order dependency structure (nei = 4).
#  # Beta <- c(0.6, 0.6, 0.6, 0.6) # Anisotropic configuration for a second
#                                  # order dependency structure (nei = 8).
#  
#  # Number of colors. Automatically set to 2 if not specified.
#  K <- 2
#  # Number of neighbors. Automatically set to 4 if not specified.
#  G <- 4
#  
#  # Optional potential on sites. Automatically set to NULL if not specified
#  potential <- runif(K,-1,1)
#  # Optional borders. Automatically set to NULL if not specified
#  Top <- Bottom <- sample(0:(K-1), width, replace = TRUE)
#  Left <- Right <- sample(0:(K-1), height, replace = TRUE)
#  Corner <- sample(0:(K-1), 4, replace = TRUE)
#  
#  # Exact sampling for the default setting
#  exact.mrf(h = height, w = width, param = Beta, view = TRUE)
#  # When specifying the number of colors and neighbors
#  exact.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta,
#            view = TRUE)
#  # When specifying an optional potential on sites
#  exact.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta,
#         pot = potential, view = TRUE)
#  # When specifying possible borders. The users will omit to mention all
#  # the non-existing borders
#  exact.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta,
#            top = Top, left = Left, bottom = Bottom,
#            right = Right, corner = Corner, view = TRUE)

## ----eval = FALSE-----------------------------------------
#  # sampler.mrf
#  # Algorithm settings
#  n <- 200
#  method <- "Gibbs"
#  
#  # Dimension of the lattice
#  height <- width <- 100
#  
#  # Interaction parameter
#  Beta <- 0.6 # Isotropic configuration
#  # Beta <- c(0.6, 0.6) # Anisotropic configuration for a first
#                        # order dependency structure (nei = 4).
#  # Beta <- c(0.6, 0.6, 0.6, 0.6) # Anisotropic configuration for a second
#                                  # order dependency structure (nei = 8).
#  
#  # Number of colors. Automatically set to 2 if not specified.
#  K <- 2
#  # Number of neighbors. Automatically set to 4 if not specified.
#  G <- 4
#  
#  # Optional potential on sites. Automatically set to NULL if not specified
#  potential <- runif(K,-1,1)
#  # Optional borders. Automatically set to NULL if not specified
#  Top <- Bottom <- sample(0:(K-1), width, replace = TRUE)
#  Left <- Right <- sample(0:(K-1), height, replace = TRUE)
#  Corner <- sample(0:(K-1), 4, replace = TRUE)
#  
#  # Sampling method for the default setting
#  sampler.mrf(iter = n, sampler = method, h = height, w = width,
#              param = Beta, view = TRUE)
#  # Sampling using an existing configuration as starting point
#  sampler.mrf(iter = n, sampler = method, h = height, w = width,
#              ncolors = K, nei = G, param = Beta,
#              initialise = FALSE, view = TRUE)
#  # Specifying optional arguments. The users may omit to mention all
#  # the non-existing borders
#  sampler.mrf(iter = n, sampler = method, h = height, w = width,
#              ncolors = K, nei = G, param = Beta,
#              pot = potential, top = Top, left = Left, bottom = Bottom,
#              right = Right, corner = Corner, view = TRUE)
#  # Gibbs sampler with sequential updates of the sites.
#  sampler.mrf(iter = n, sampler = "Gibbs", h = height, w = width,
#              ncolors = K, nei = G, param = Beta,
#              random = FALSE, view = TRUE)

