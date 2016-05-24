# Copyright (C) 2015 Johannes Ranke
# Contact: jranke@uni-bremen.de
# The summary function is an adapted and extended version of summary.modFit
# from the FME package, v 1.1 by Soetart and Petzoldt, which was in turn
# inspired by summary.nls.lm

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>

mmkin <- function(models = c("SFO", "FOMC", "DFOP"), datasets, 
                  cores = round(detectCores()/2), cluster = NULL, ...) 
{
  parent_models_available = c("SFO", "FOMC", "DFOP", "HS", "SFORB", "IORE") 
  n.m <- length(models)
  n.d <- length(datasets)
  n.fits <- n.m * n.d
  fit_indices <- matrix(1:n.fits, ncol = n.d)

  # Check models and define their names
  if (!all(sapply(models, function(x) inherits(x, "mkinmod")))) {
    if (!all(models %in% parent_models_available)) {
      stop("Please supply models as a list of mkinmod objects or a vector combined of\n  ",
           paste(parent_models_available, collapse = ", ")) 
    } else {
      names(models) <- models
    } 
  } else {
    if (is.null(names(models))) names(models) <- as.character(1:n.m)
  }

  # Check datasets and define their names
  if (is.null(names(datasets))) names(datasets) <- as.character(1:n.d)

  # Define names for fit index
  dimnames(fit_indices) <- list(model = names(models),
                                dataset = names(datasets))


  fit_function <- function(fit_index) {
    w <- which(fit_indices == fit_index, arr.ind = TRUE)
    model_index <- w[1]
    dataset_index <- w[2]
    mkinfit(models[[model_index]], datasets[[dataset_index]], ...)
  }

  if (is.null(cluster)) {
    results <- mclapply(as.list(1:n.fits), fit_function, mc.cores = cores)
  } else {
    results <- parLapply(cluster, as.list(1:n.fits), fit_function)
  }

  attributes(results) <- attributes(fit_indices)
  class(results) <- "mmkin"
  return(results)
}

"[.mmkin" <- function(x, i, j, ..., drop = FALSE) {
  class(x) <- NULL
  x_sub <- x[i, j, drop = drop]
  if (!drop) class(x_sub) <- "mmkin"
  return(x_sub)
}
