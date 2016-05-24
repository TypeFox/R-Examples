# Author: Richard Bourgon (coauthor of package 'genefilter')
# Algorithm is due to D.R. Bickel
hrm <-
function(x,         # sample (the data)
         bw = NULL, # bandwidth (fraction of the observations to consider)
         ...)
         #! introduire l'argument 'k'
{
############################################################
# Bickel's mode estimator
# HRM = half-range mode
# C routine 'half_range_mode' comes from package 'genefilter'
############################################################

  if (is.null(bw)) stop("argument 'bw' is missing")
  if (bw <= 0 | bw > 1) stop("argument 'bw' must belong to (0, 1]")
  
  y <- sort(x)    
  return(.C("half_range_mode",
            data = as.double(y),
            n = as.integer(length(y)),
            beta = as.double(bw),
            diag = as.integer(FALSE),
            M = double(1),
            PACKAGE = "modeest")$M)
} 

#half.range.mode <- Bickel <- bickel <- HRM <- hrm
