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

hergm.mcmc <- function(original.formula, max_number, initialize, network, model, hyper_prior, parametric, MHproposal, MCMCparams, verbose, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, mean_between, eta, indicator, parallel, simulate, seeds, mh_scale, variational, temperature, predictions, perturb) 
{
  # Prepare
  if (simulate == FALSE) scalefactor <- hergm.set.mcmc(max_number, initialize, network, model, hyper_prior, parametric, MHproposal, MCMCparams, verbose, indicator, alpha_shape, alpha_rate, alpha,  eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, mean_between, eta, simulate, parallel, seeds, predictions, variational, temperature, mh_scale, perturb) # The last argument is the initial value of the scale factor
  Clist <- ergm.Cprepare(network, model)
  if (Clist$dir == FALSE) maxedges <- Clist$n * (Clist$n - 1) / 2 # Undirected
  else maxedges <- Clist$n * (Clist$n - 1) # Directed
  hergm_list <- hergm.preprocess(max_number, initialize=FALSE, network, model, hyper_prior, parametric, Clist, MHproposal, MCMCparams, maxedges, alpha_shape, alpha_rate, alpha,  eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, mean_between, eta, indicator, simulate, parallel, variational, temperature, predictions, verbose, perturb)
  if (simulate == FALSE) hergm_list$scalefactor <- scalefactor # Set scale factor
  else if (hergm_list$hyper_prior == 1)
    {
    k <- number_clusters(hergm_list$alpha, hergm_list$Clist$n)
    cat("\nDirichlet process:")
    cat("\n- scaling parameter: ", hergm_list$alpha, sep = "")
    cat("\n- mean of number of non-empty blocks: ", formatC(k$mean, digits = 2, width = 4, format = "f", mode = "real"), sep = "") 
    cat("\n- variance of number of non-empty blocks: ", formatC(k$variance, digits = 2, width = 4, format = "f", mode = "real"), sep = "") 
    cat("\n")   
    }
  flush.console() # Windows

  #print(model$term)
  #ergm.getglobalstats(network,model)
  #model<-append.rhs.formula(model,list(as.name("edges")))
  #model<-append.rhs.formula(model,list(as.name("edges")))
  #print(model$term)
  #hergm_list$Clist <- ergm.Cprepare(network, model)
  #print("done")

  # Run
  if (simulate == TRUE) parallel <- 1
  if (parallel > 1) # Parallel computing
    {
    number <- parallel # Specify number of computing nodes on cluster
    if (length(seeds) < number) 
      {
      cluster.seeds <- runif(number)
      }
    else 
      {
      cluster.seeds <- seeds
      }
    cluster <- makePSOCKcluster(number)
    clusterEvalQ(cluster, library(hergm))
    s <- clusterApplyLB(cluster, cluster.seeds[1:number], hergm.wrapper, hergm_list)
    stopCluster(cluster)
    mcmc <- append(s[[1]]$mcmc, s[[2]]$mcmc)
    if (number > 2) 
      {
      for (i in 3:number) mcmc <- append(mcmc, s[[i]]$mcmc)
      }
    sample <- list()
    sample$mcmc <- mcmc
    }
  else sample <- hergm.wrapper(seeds[1], hergm_list) # Non-parallel computing

  # Store
  output <- list() 
  output$network <- network
  output$n <- Clist$n
  output$model <- model
  output$max_number <- hergm_list$max_number
  output$number_fixed <- hergm_list$number_fixed
  output$d1 <- hergm_list$d1
  output$d2 <- hergm_list$d2
  output$parallel <- hergm_list$parallel
  output$simulate <- hergm_list$simulate
  output$sample_size <- min(10000, hergm_list$MCMCparams$samplesize)
  if (simulate == TRUE) 
    {
    output$sample <- sample$sample
    number_edges <- sample$sample_heads[1]
    if (number_edges == 0)
      {
      sample_heads <- NULL
      sample_tails <- NULL
      }
    else
      { 
      sample$sample_heads <- sample$sample_heads[-1]
      sample$sample_tails <- sample$sample_tails[-1]
      output$heads <- sample$sample_heads[1:number_edges]
      output$tails <- sample$sample_tails[1:number_edges]
      }
    }
  output$predictions <- predictions
  output$sample <- sample$mcmc 

  output
}

