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

hergm.set.mcmc <- function(max_number, initialize, network, model, hyper_prior, parametric, MHproposal, MCMCparams, verbose, indicator, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, mean_between, simulate, parallel, seeds, predictions, variational, temperature, mh_scale, perturb)
{

  # verbose <- 1

  # Prepare I
  cp_samplesize <- MCMCparams$samplesize # Store
  MCMCparams$samplesize <- round(parallel * cp_samplesize / 100)
  if (MCMCparams$samplesize < 100) MCMCparams$samplesize <- 100
  else if (MCMCparams$samplesize > 1000) MCMCparams$samplesize <- 1000

  # Prepare II
  Clist <- ergm.Cprepare(network, model)
  if (Clist$dir == FALSE) maxedges <- Clist$n * (Clist$n - 1) / 2 # Undirected
  else maxedges <- Clist$n * (Clist$n - 1) # Directed
  hergm_list <- hergm.preprocess(max_number, initialize, network, model, hyper_prior, parametric, Clist, MHproposal, MCMCparams, maxedges, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, mean_between, eta, indicator, simulate, parallel = 1, variational, temperature, predictions = FALSE, verbose = -1, perturb)
  # Metropolis-Hastings: finding scale factor
  if (verbose >= 0) 
    {
    #if (hergm_list$dependence > 0) cat("\nMCMC: mean-field methods generate candidates of block memberships.\n")
    cat("\nMetropolis-Hastings algorithm:")
    cat("\n   ", formatC("ergm", digits = 0, width = 16, format = "f", mode = "character"), sep = "")
    }
  if (hergm_list$dependence == 0) 
    {
    number <- 2
    if (verbose >= 0) cat(formatC("hergm", digits = 0, width = 8, format = "f", mode = "character"), sep = "")
    }
  else 
    {
    number <- 2 + (Clist$n - hergm_list$min_size)
    if (verbose >= 0)
      {
      for (i in hergm_list$min_size:Clist$n)
        {
        cat(formatC(i, digits = 0, width = 8, format = "f", mode = "integer"), sep = "")
        }   
      }
    }
  if (is.null(mh_scale) || (mh_scale < 0)) scalefactor <- rep.int(1, number)
  else scalefactor <- rep.int(mh_scale, number)
  hergm_list$scalefactor <- scalefactor
  min_accept <- 0.2
  iteration <- 1
  s <- hergm.wrapper(seeds[1], hergm_list)
  if (verbose >= 0) 
    {
    cat("\n", "(", formatC(iteration, digits = 0, width = 1, mode = "character"), ")", sep = "")
    cat(formatC("scale:", digits = 0, width = 8, format = "f", mode = "character"), sep = "")
    for (i in 1:number) 
      {
      cat(formatC(scalefactor[i], digits = 4, width = 8, format = "f", mode = "real"), sep = "")
      }
    cat("\n   ", formatC("rate:", digits = 0, width = 8, format = "f", mode = "character"), sep = "")
    for (i in 1:number) 
      {
      cat(formatC(s$mh_accept[i], digits = 4, width = 8, format = "f", mode = "real"), sep = "")
      }
    }
  while ((s$mh_accept[1] < min_accept) && (min(s$mh_accept[2:number]) < (min_accept / 2)) && (iteration <= 20))
    { 
    iteration <- iteration + 1
    hergm_list <- hergm.preprocess(number, initialize, network, model, hyper_prior, parametric, Clist, MHproposal, MCMCparams, maxedges, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, mean_between, eta, indicator, simulate, parallel = 1, variational, temperature, predictions = FALSE, verbose = -1, perturb)
    for (i in 1:number) 
      {
      if (s$mh_accept[i] < min_accept) scalefactor[i] <- scalefactor[i] / 2
      if (i > 2) 
        {
        if (scalefactor[i] > scalefactor[i-1]) scalefactor[i] <- scalefactor[i-1]
        }
      }
    hergm_list$scalefactor <- scalefactor
    s <- hergm.wrapper(seeds[1], hergm_list)
    if (verbose >= 0) 
      {
      cat("\n", "(", formatC(iteration, digits = 0, width = 1, mode = "character"), ")", sep = "")
      cat(formatC("scale:", digits = 0, width = 8, format = "f", mode = "character"), sep = "")
      for (i in 1:number) 
        {
        cat(formatC(scalefactor[i], digits = 4, width = 8, format = "f", mode = "real"), sep = "")
        }
      cat("\n   ", formatC("rate:", digits = 0, width = 8, format = "f", mode = "character"), sep = "")
      for (i in 1:number) 
        {
        cat(formatC(s$mh_accept[i], digits = 4, width = 8, format = "f", mode = "real"), sep = "")
        }
      }
    }

  if (verbose >= 0) cat("\n")

  MCMCparams$samplesize <- cp_samplesize # Reset

  scalefactor

}

