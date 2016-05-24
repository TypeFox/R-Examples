# Copyright (C) 2015
# Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.



# --------------------------------------
# rgcca for multiple integration  and a regularisation parameter tau
#--------------------------------------
wrapper.rgcca <- function (blocks, 
                           design =  1 - diag(length(blocks)), 
                           ncomp = rep(2, length(blocks)), 
                           tau = "optimal",
                           scheme = "centroid", 
                           scale = TRUE, 
                           bias = FALSE, 
                           max.iter = 1000,
                           tol = .Machine$double.eps, 
                           verbose = FALSE, 
                           near.zero.var = FALSE
) {
  
  # note here: few hard coded arguments:
  # mode  = canonical, penalty = NULL, keep = NULL, indY = NULL, init (Mode of initialization use in the SGCCA algorithm with Singular Value Decompostion ("svd") 
  result.rgcca = srgcca(blocks = blocks, indY = NULL, design = design, tau = tau, 
                        ncomp = ncomp, scheme = scheme, scale = scale, bias = bias, init = "svd.rgcca", 
                        tol = tol, verbose = verbose, mode = "canonical", max.iter = max.iter, 
                        keep = NULL, near.zero.var = near.zero.var, penalty = NULL)
  
  result.rgcca$class = "rgcca"
  class(result.rgcca) = "rgcca"
  return(invisible(result.rgcca))
}


# --------------------------------------
# sgcca for multiple integration with variable selection in each block (lasso penalisation)
#--------------------------------------
wrapper.sgcca <- function (blocks, 
                           design =  1 - diag(length(blocks)),
                           penalty = NULL, 
                           ncomp = rep(2, length(blocks)),
                           keep = NULL, 
                           scheme = "centroid",
                           scale = TRUE, 
                           bias = FALSE, 
                           max.iter = 1000,
                           tol = .Machine$double.eps, 
                           verbose = FALSE, 
                           near.zero.var = FALSE
                           ) {
  
  
    # note here: mode is hard coded as canonical
  result.sgcca = srgcca(blocks = blocks, indY = NULL, design = design, tau = rep(1, length(blocks)), 
                        ncomp = ncomp, scheme = scheme, scale = scale, bias = bias, init = "svd.sgcca", 
                        tol = tol, verbose = verbose, mode = "canonical", max.iter = max.iter, 
                        keep = keep, near.zero.var = near.zero.var, penalty = penalty)

  result.sgcca$class = "sgcca"
  class(result.sgcca) = "sgcca"
  return(invisible(result.sgcca))
}



# --------------------------------------
# srgccda for multiple integration with a block (or more) as outcome(s) in Y and variable selection (lasso penalisation)
#--------------------------------------
wrapper.sgccda <- function(
                            blocks, 
                            Y, 
                            design = NULL, 
                            ncomp = rep(2, length(blocks)) ,                           
                            keep = NULL, 
                            scheme = "centroid", 
                            scale = TRUE, 
                            bias = FALSE, 
                            max.iter = 1000,
                           tol = .Machine$double.eps, 
                           verbose = FALSE, 
                           near.zero.var = FALSE
                           ) {

  #-- Define blocks
  if (!is.list(blocks))
    stop("'blocks' must be a list containing the data sets.", call. = FALSE)  
  
  if (is.null(dim(Y))) {
    Y = as.factor(Y)
    ind.mat = data.frame(unmap(as.numeric(Y)))
    colnames(ind.mat) = levels(Y)
  } else {
    stop("'Y' should be a factor or a class vector.")
  }
  
  #-- ncomp
  if (!is.vector(ncomp, mode = "numeric"))
    stop("'ncomp' must be a vector")
  
  if (length(ncomp) == length(blocks)){
    ncomp = c(ncomp, max(ncomp))
    message("'ncomp' has changed to include Y as a block")
  }
  
  #-- Design
  if (is.null(design)) {
    design = 1 - diag(length(blocks) + 1)
  } else if (ncol(design) != nrow(design) || ncol(design) < length(blocks) || ncol(design) > (length(blocks) + 1) || any(!design %in% c(0,1))){
    stop('Invalid design matrix.')
  } else if (ncol(design) == length(blocks)){ 
    message('Design matrix has changed to include Y as a block')
    design = rbind(cbind(design, 1), 1)
    diag(design) = 0
  }
  
  #-- keep
  if (is.list(keep) & (length(keep) == length(blocks))){
    keep[[length(keep) + 1]] = rep(nlevels(Y), ncomp[length(blocks) + 1])
    message("'keep' has changed and include Y as a block")
  }
  
  # merge blocks and Y only to run srgcca, the output "blocks" will not contain Y
  blocks[[length(blocks) + 1]] = ind.mat
  names(blocks)[length(blocks)] = "Y"
    
  # note here: mode is hard coded as regression as we are performing a supervised analysis w.r.t Y
  result.sgccada = srgcca(blocks = blocks, indY = length(blocks), design = design, tau = rep(1, length(blocks)),
                        ncomp = ncomp, scheme = scheme, scale = scale, bias = bias, init = "svd.da", 
                        tol = tol, verbose = verbose, mode = "regression", max.iter = max.iter,
                        keep = keep, near.zero.var = near.zero.var, penalty = NULL)
  
  result.sgccada$Y = Y
  result.sgccada$ind.mat = ind.mat
  result.sgccada$class = c("sgccda","sgcca")
  class(result.sgccada) = c("sgccda","sgcca")
  return(invisible(result.sgccada))
}