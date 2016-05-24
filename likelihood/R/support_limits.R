########################################################
# support_limits
# This function calculates support limits for those parameters that were
# optimized through the simulated annealing process.  The algorithm is
# just to walk out the function until the limits are reached.  This function
# will only work perfectly if the likelihood function has a smooth surface
# and can be stymied by bumps.  However, it's still pretty darn useful.
#
# Arguments:
# model - the function which contains the model statement whose parameters
# were optimized.
# par - the parameter list as it was used in likeli and/or neighanneal,
# with the best values for the varying parameters.
# var - list of other variables and data needed by the model
# and pdf functions, any type as needed.
# source_data - the source data as it was passed to likeli and/or
# neighanneal.
# pdf - the probability distribution function as it was passed to likeli
# and/or neighanneal.
# par_lo - list containing the lower limits of the varying parameters.  Can
# be a ragged array but no more than 2 dimensions.  The names of the list
# members must match those in par. Any values ommitted are assumed to be
# negative infinity.
# par_hi - list containing the upper limits of the varying parameters.  Can
# be a ragged array but no more than 2 dimensions.  The names of the list
# members must match those in par. Any values ommitted are assumed to be
# infinity.
# delta - the number of pieces into which to divide each parameter.  This
# is the size of the "step" the function takes in trying to find the support
# limits.
# slimit - the number of likelihood units less than the optimum
# likelihood for which the support intervals will be calculated.  So if
# the max likelihood is -98, then for each parameter we change it in steps
# of delta until we reach a likelihood of -100.
#
# This is based on code originally written in Delphi by Charles D. Canham.
# Author:  Lora Murphy, Cary Institute of Ecosystem Studies
# murphyl@caryinstitute.org
########################################################
support_limits<-function(model, par, var, source_data, pdf, par_lo = NULL,
   par_hi = NULL, delta = 100, slimit = 2) {

  # Check that par_lo and par_hi have corresponding labels in par,
  # and replace missing values with "infinity"
  if (is.null(par_lo)) {
    par_lo <- par
    for (i in 1:length(par)) {
      par_lo[[i]] = rep(-.Machine$double.xmax, times = length(par[[i]]))
    }
  }
  if (is.null(par_hi)) {
    par_hi <- par
    for (i in 1:length(par)) {
      par_hi[[i]] = rep(.Machine$double.xmax, times = length(par[[i]]))
    }
  }
  for (i in 1:length(par)) {
    if (!any(names(par_lo) == names(par)[i])) {
      par_lo[[names(par)[i]]] = rep(-.Machine$double.xmax, times = length(par[[i]]))
    }
    if (!any(names(par_hi) == names(par)[i])) {
      par_hi[[names(par)[i]]] = rep(.Machine$double.xmax, times = length(par[[i]]))
    }
  }

  # Replace any true values of infinity with our fake infinity
  for (i in 1:length(par_lo)) {
    par_lo[[i]] = ifelse(is.infinite(par_lo[[i]]), -.Machine$double.xmax, par_lo[[i]])
  }
  for (i in 1:length(par_hi)) {
    par_hi[[i]] = ifelse(is.infinite(par_hi[[i]]), .Machine$double.xmax, par_hi[[i]])
  }

  lower_limit<-par # where we'll stash lower limits
  upper_limit<-par # where we'll stash lower limits

  # Get the best likelihood
  best_lh <- likeli(model, par, var, source_data, pdf)
  if (is.infinite(best_lh) || is.nan(best_lh) || is.na(best_lh)) {
    for (i in 1:length(lower_limit)) lower_limit[[i]] <- NA
    for (i in 1:length(upper_limit)) upper_limit[[i]] <- NA
    return(list(upper_limits = upper_limit, lower_limits = lower_limit))
   #stop("Cannot calculate support limits:  best likelihood is infinite or not a number.\n")
  }

  # Get a fake likelihood value for replacing infinite or NaN values
  # that will be beyond the support limits
  bad_likeli <-  best_lh - (slimit + 10)

  # Loop through the varying parameters - check for the possibility of
  # ragged arrays
  for (i in 1:length(par)) {
    for (j in 1:length(par[[i]])) {

      # Make a copy of the parameter we're about to vary so we can set
      # it back to the optimum at the end
      par_copy <- par[[i]][[j]]

      smallest_step <- par[[i]][[j]]/delta

      ################################
      ## Calculate the upper limit
      ################################

      # First - check likelihood at upper limit - if it's within the support
      # interval we don't have to do anything else
      par[[i]][[j]] <- par_hi[[i]][[j]]
      lhood <- likeli(model, par, var, source_data, pdf)
      if (is.infinite(lhood) || is.nan(lhood) || is.na(lhood)) lhood <- bad_likeli
      lhdiff <- best_lh - lhood

      if (!is.nan(lhdiff) && !is.na(lhdiff) && lhdiff < slimit) {
        upper_limit[[i]][[j]] <- par_hi[[i]][[j]]
      }
      else {
        # Set the parameter back to the optimum
        par[[i]][[j]] <- par_copy
        lhdiff <- 0

        # Set up the biggest step size - 1/50th the search range
        biggest_step <- (par_hi[[i]][[j]] - par[[i]][[j]])/50
        if(is.infinite(biggest_step)) biggest_step <- 1e306

        # Perform the search until we've refined the search as much as possible
        step <- biggest_step
        while (step > smallest_step) {

          # Step along at the current step size
          lhdiff<-0
          while (!is.nan(lhdiff) && !is.na(lhdiff) && (lhdiff < slimit) && (par[[i]][[j]] < par_hi[[i]][[j]])) {
            par[[i]][[j]] <- par[[i]][[j]] + step
            lhood <- likeli(model, par, var, source_data, pdf)
            if (is.infinite(lhood) || is.nan(lhood) || is.na(lhood)) lhood <- bad_likeli
            lhdiff <- best_lh - lhood
          }

          # Step back once so we can search more finely on
          # the appropriate interval on the next loop
          par[[i]][[j]] <- par[[i]][[j]] - step
          # Reduce the step size by one order of magnitude
          step <- step / 10
          # It is possible to have a step size so small that adding it produces
          # a number indistinguishable from the previous value. If this
          # is the case, quit.
          if (!((par[[i]][[j]] + step) > par[[i]][[j]])) {
            break
          }
        }

        # Do the last search on the smallest step size, if it's not already
        # too small
        step <- smallest_step
        if ((par[[i]][[j]] + step) > par[[i]][[j]]) {
          lhdiff<-0
          while (!is.nan(lhdiff) && !is.na(lhdiff) && (lhdiff < slimit) && (par[[i]][[j]] < par_hi[[i]][[j]])) {
            par[[i]][[j]] <- par[[i]][[j]] + step
            lhood <- likeli(model, par, var, source_data, pdf)
            if (is.infinite(lhood) || is.nan(lhood) || is.na(lhood)) lhood <- bad_likeli
            lhdiff <- best_lh - lhood
          }
        }
        
        upper_limit[[i]][[j]] <- par[[i]][[j]] - step

        # Set the parameter back to the optimum
        par[[i]][[j]] <- par_copy
      }

      ################################
      ## Calculate the lower limit
      ################################

      # First - check likelihood at lower limit - if it's within the support
      # interval we don't have to do anything else
      par[[i]][[j]] <- par_lo[[i]][[j]]
      lhood <- likeli(model, par, var, source_data, pdf)
      if (is.infinite(lhood) || is.nan(lhood) || is.na(lhood)) lhood <- bad_likeli
      lhdiff <- best_lh - lhood

      if (!is.nan(lhdiff) && !is.na(lhdiff) && lhdiff < slimit) {
        lower_limit[[i]][[j]] <- par_lo[[i]][[j]]
      }
      else {
        # Set the parameter back to the optimum
        par[[i]][[j]] <- par_copy
        lhdiff <- 0

        # Set up the biggest step size - 1/50th the search range
        biggest_step <- (par[[i]][[j]] - par_lo[[i]][[j]])/50

        # Perform the search until we've refined the search as much as possible
        step <- biggest_step
        if(is.infinite(biggest_step)) biggest_step <- 1e306
        while (step > smallest_step) {

          # Step along at the current step size
          lhdiff<-0
          while (!is.nan(lhdiff) && !is.na(lhdiff) && (lhdiff < slimit) && (par_lo[[i]][[j]] < par[[i]][[j]])) {
            par[[i]][[j]] <- par[[i]][[j]] - step
            lhood <- likeli(model, par, var, source_data, pdf)
            if (is.infinite(lhood) || is.nan(lhood) || is.na(lhood)) lhood <- bad_likeli
            lhdiff <- best_lh - lhood
          }

          # Step back once so we can search more finely on
          # the appropriate interval on the next loop
          par[[i]][[j]] <- par[[i]][[j]] + step
          # Reduce the step size by one order of magnitude
          step <- step / 10
          # It is possible to have a step size so small that subtracting it 
          # produces a number indistinguishable from the previous value. If
          # this is the case, quit.
          if (!((par[[i]][[j]] - step) < par[[i]][[j]])) {
            break
          }
        }

        # Do the last search on the smallest step size
        step <- smallest_step
        if ((par[[i]][[j]] - step) < par[[i]][[j]]) {
          lhdiff<-0
          while (!is.nan(lhdiff) && !is.na(lhdiff) && (lhdiff < slimit) && (par_lo[[i]][[j]] < par[[i]][[j]])) {
            par[[i]][[j]] <- par[[i]][[j]] - step
            lhood <- likeli(model, par, var, source_data, pdf)
            if (is.infinite(lhood) || is.nan(lhood) || is.na(lhood)) lhood <- bad_likeli
            lhdiff <- best_lh - lhood
          }
        }
          
        lower_limit[[i]][[j]] <- par[[i]][[j]]+step

        # Set the parameter back to the optimum
        par[[i]][[j]] <- par_copy
      }
    } # end of parameter loop
  }

  list(upper_limits = upper_limit, lower_limits = lower_limit)
}

