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

hergm.wrapper <- function(seed, hergm_list) 
{
  # print(hergm_list)
  if (is.null(seed) == FALSE) set.seed(seed)
  if (hergm_list$simulate == TRUE)
    {
    sample <- .C("Simulation",
    as.integer(hergm_list$dependence),
    as.integer(hergm_list$hierarchical),
    as.integer(hergm_list$d),
    as.integer(hergm_list$d1),
    as.integer(hergm_list$d2),
    as.integer(hergm_list$structural),
    as.integer(hergm_list$min_size),
    as.integer(hergm_list$max_number),
    as.integer(hergm_list$null$alpha),
    as.integer(hergm_list$null$eta),
    as.integer(hergm_list$null$indicator),
    as.double(hergm_list$alpha),
    as.double(hergm_list$alpha_shape),
    as.double(hergm_list$alpha_rate),
    as.double(hergm_list$eta_mean1),
    as.double(hergm_list$eta_mean2),
    as.double(hergm_list$cf1),
    as.double(hergm_list$cf2),
    as.double(hergm_list$p1),
    as.double(hergm_list$p2),
    as.double(hergm_list$eta_mean_mean),
    as.double(hergm_list$eta_mean_precision),
    as.double(hergm_list$eta_precision_shape),
    as.double(hergm_list$eta_precision_rate),
    as.double(hergm_list$eta),
    as.integer(hergm_list$indicator),
    as.integer(hergm_list$Clist$heads), 
    as.integer(hergm_list$Clist$tails),
    as.integer(hergm_list$Clist$nedges), 
    as.integer(hergm_list$MCMCparams$MCMC.init.maxedges), 
    as.integer(hergm_list$Clist$n),
    as.integer(hergm_list$Clist$dir), 
    as.integer(hergm_list$Clist$bipartite),
    as.integer(hergm_list$Clist$nterms), 
    as.character(hergm_list$Clist$fnamestring), 
    as.character(hergm_list$Clist$snamestring),
    as.character(hergm_list$MHproposal$name), 
    as.character(hergm_list$MHproposal$pkgname),
    as.double(hergm_list$Clist$inputs),
    as.double(hergm_list$Clist$inputs),
    as.double(hergm_list$theta),
    as.integer(hergm_list$MCMCparams$samplesize),
    sample = as.double(t(hergm_list$MCMCparams$stats)),
    as.integer(hergm_list$MCMCparams$burnin), 
    as.integer(hergm_list$MCMCparams$interval),
    as.integer(hergm_list$verbose),
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$attribs),
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$maxout), 
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$maxin),
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$minout), 
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$minin),
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$condAllDegExact), 
    as.integer(length(hergm_list$MHproposal$arguments$constraints$bd$attribs)),
    as.integer(hergm_list$maxedges),
    as.integer(hergm_list$max_iteration),
    as.integer(hergm_list$between),
    as.double(hergm_list$mean_between),
    as.integer(hergm_list$predictions),
    mcmc = as.double(hergm_list$mcmc),
    sample_heads = as.integer(hergm_list$sample_heads),
    sample_tails = as.integer(hergm_list$sample_tails),
    as.integer(hergm_list$prior_assumptions),
    status = integer(1), 
    PACKAGE="hergm")
    }
  else
    {
    sample <- .C("Inference",
    as.integer(hergm_list$model_type),
    as.integer(hergm_list$dependence),
    as.integer(hergm_list$hierarchical),
    as.integer(hergm_list$decomposable),
    as.integer(hergm_list$d),
    as.integer(hergm_list$d1),
    as.integer(hergm_list$d2),
    as.integer(hergm_list$structural),
    as.integer(hergm_list$min_size),
    as.integer(hergm_list$max_number),
    as.double(hergm_list$alpha),
    as.double(hergm_list$alpha_shape),
    as.double(hergm_list$alpha_rate),
    as.double(hergm_list$eta_mean1),
    as.double(hergm_list$eta_mean2),
    as.double(hergm_list$cf1),
    as.double(hergm_list$cf2),
    as.double(hergm_list$p1),
    as.double(hergm_list$p2),
    as.double(hergm_list$eta_mean_mean),
    as.double(hergm_list$eta_mean_precision),
    as.double(hergm_list$eta_precision_shape),
    as.double(hergm_list$eta_precision_rate),
    as.integer(hergm_list$indicator),
    as.integer(hergm_list$Clist$heads), 
    as.integer(hergm_list$Clist$tails),
    as.integer(hergm_list$Clist$nedges), 
    as.integer(hergm_list$MCMCparams$MCMC.init.maxedges), 
    as.integer(hergm_list$Clist$n),
    as.integer(hergm_list$Clist$dir), 
    as.integer(hergm_list$Clist$bipartite),
    as.integer(hergm_list$Clist$nterms), 
    as.character(hergm_list$Clist$fnamestring), 
    as.character(hergm_list$Clist$snamestring),
    as.character(hergm_list$MHproposal$name), 
    as.character(hergm_list$MHproposal$pkgname),
    as.double(hergm_list$Clist$inputs),
    as.double(hergm_list$Clist$inputs),
    as.double(hergm_list$theta),
    as.integer(hergm_list$MCMCparams$samplesize),
    sample = as.double(t(hergm_list$MCMCparams$stats)),
    as.integer(hergm_list$MCMCparams$burnin), 
    as.integer(hergm_list$MCMCparams$interval),
    newnwheads = integer(hergm_list$maxedges),
    newnwtails = integer(hergm_list$maxedges),
    as.integer(hergm_list$verbose), 
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$attribs),
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$maxout), 
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$maxin),
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$minout), 
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$minin),
    as.integer(hergm_list$MHproposal$arguments$constraints$bd$condAllDegExact), 
    as.integer(length(hergm_list$MHproposal$arguments$constraints$bd$attribs)),
    as.integer(hergm_list$maxedges),
    as.integer(hergm_list$max_iteration),
    as.integer(hergm_list$between),
    as.integer(hergm_list$predictions),
    mcmc = as.double(hergm_list$mcmc),
    as.double(hergm_list$scalefactor),
    mh_accept = as.double(hergm_list$mh_accept),
    as.double(hergm_list$q_i),
    as.integer(hergm_list$parallel),
    as.double(hergm_list$temperature),
    as.integer(hergm_list$prior_assumptions),
    status = integer(1),
    PACKAGE="hergm")
    }
  sample
}

