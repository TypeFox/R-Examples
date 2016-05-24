# Check the given X and Y limit list. Each list can have 2 elements:
# XSpace and YSpace. If the element is NULL, it is set, otherwise it is
# checked. If the dimensionality of X and Y space is greater than 2, the limits
# are set to NULL. The same happens if there is a discrete variable in the 2D case.
# We don't need limits in this cases. 
# If the dimensionality is 1 for X or Y Space, for this space the ylim is NULL.
getOptPathLims = function(xlim, ylim, op.x, op.y, iters, scale) {
  
  assertList(xlim)
  assertList(ylim)
  
  # First, we calculate for XSpace and YSpace plot the limits
  # Nasty ifs, since this can be one of half a dozen plots, and for every
  # plot we have different defaults
  for (space in c("XSpace", "YSpace")) {
    op.frame = if (space == "XSpace") op.x else op.y
    dim = ncol(op.frame)
    classes = BBmisc::vcapply(op.frame, function(x) class(x))
    
    # For Multi-D Plot, no limits are needed. Warn, if the user specified some
    # and set to NULL
    if (dim > 2L) {
      if (!is.null(xlim[[space]]))
        warning(paste("You specified xlims for multi-D plot in",
          space, "but xlims for this plots are not supported."))
      if (!is.null(ylim[[space]]))
        warning("You specified ylims for multi-D plot in",
          space, "but ylims for this plots are not supported.")
      xlim[[space]] = NULL
      ylim[[space]] = NULL
      next
    }
    
    # For 1Dnumeric this is easy - either check user input
    # or set to min - scale * range and max + scale * range 
    if (dim == 1L && (classes == "numeric")) {
      if (is.null(xlim[[space]])) {
        xlim[[space]] = range(op.frame[, 1L])
        xlim[[space]] = c(-1, 1) * scale * abs(diff(xlim[[space]])) + xlim[[space]]
      } else {
        assertNumeric(xlim[[space]], len = 2L, any.missing = FALSE)
      }
      if (!is.null(ylim[[space]]))
        warning("You specified ylims for 1D numeric plot in",
          space, "but ylims for this plots are not supported.")
      ylim[[space]] = NULL
      next
    }
    
    # limits for barplot (1D discrete case)
    # Here, xlims are not meaningful, ylim is the modal value
    if (dim == 1L && classes == "factor") {
      if (!is.null(xlim[[space]]))
        warning(paste("You specified xlims for 1D barplot in",
          space, "but xlims for this plots are not supported."))
      xlim[[space]] = NULL
      
      if (is.null(ylim[[space]])) {
        ylim[[space]] = c(0, max(table(op.frame)))
      } else {
        assertNumeric(ylim[[space]], len = 2L, any.missing = FALSE)
      }
      next
    }
    
    # For dim = 2L we have to check for both variables, if they are discrete
    # for discrete, we don't want limits, for factors they are not meaningful
    if (dim == 2L) {
      if (classes[1L] == "numeric") {
        if (is.null(xlim[[space]])) {
          xlim[[space]] = range(op.frame[, 1L])
          xlim[[space]] = c(-1, 1) * scale * abs(diff(xlim[[space]])) + xlim[[space]]
        } else {
          assertNumeric(xlim[[space]], len = 2L, any.missing = FALSE)
        }
      } else {
        if (!is.null(xlim[[space]]))
          warning(paste("You specified xlims for 2D scatter plot in",
            space, "but the variable is discrete here and therefor xlims are not supported."))
        xlim[[space]] = NULL
      }
      
      if (classes[2L] == "numeric") {
        if (is.null(ylim[[space]])) {
          ylim[[space]] = range(op.frame[, 2L])
          ylim[[space]] = c(-1, 1) * scale * abs(diff(ylim[[space]])) + ylim[[space]]
        } else {
          assertNumeric(ylim[[space]], len = 2L, any.missing = FALSE)
        }
      } else {
        if (!is.null(ylim[[space]]))
          warning(paste("You specified ylims for 2D scatter plot in",
            space, "but the variable is discrete here and therefor ylims are not supported."))
        ylim[[space]] = NULL
      } 
    }
    
  }
  
  # For now I say: The code for checking and limits for the over.time.plots
  # is way to bad and very complicated. It's not worth the afford atm, since in
  # 99% the defaults of ggplot are the way to go.
  # the user can get the plots via render, he has to set the limits here
  # for himself
  # NOTE: This code is not finished!
  
#   # Now we need the limits for over.time plots. This is a bit complicated,
#   # since we can have a list of over.time plots, so we also can get
#   # a list of lims.
#   for (space in c("x.over.time", "y.over.time")) {
#     
#     over.time.vars = 
#     
#     # First, ensure we have lists.
#     if (!is.list(xlim[[space]]))
#       xlim[[space]] = list(xlim[[space]])
#     
#     if (!is.list(ylim[[space]]))
#       ylim[[space]] = list(ylim[[space]])
#     
#     # Here, you allways have to specify limits for every plot - if you specify
#     # some limits.
#     assertList(ylim[[space]], len = )
#     assertList(ylim[[space]], 
#     
#     for (i in seq_along(get(space))) {
#       
#       op.frame = if (space == "x.over.time") op.x else op.y
#       var.names = if (space == "x.over.time") x.over.time[[i]] else y.over.time[[i]]
#       op.frame = op.frame[, var.names, drop = FALSE]
#       
#       # lim.x is iteration number here. If NULL, use ggplot defaults, else check.
#       if (is.null(xlim[[space]][[i]])) {
#         xlim[[space]][[i]] = c(NA_real_, NA_real_)
#       } else {
#         asInteger(xlim[[space]][[i]], len = 2)
#       }
#       
#       # lim.y as minimum and maximum of all plotted variables:
#       if (is.null(ylim[[space]][[i]])) {
#         ylim[[space]][[i]] = range(op.frame)
#         ylim[[space]][[i]] = c(-1, 1) * scale * abs(diff(ylim[[space]][[i]])) + ylim[[space]][[i]]
#       } else {
#         assertNumeric(ylim[[space]][[i]], len = 2L, any.missing = FALSE)
#       }
#       
#     }
#   }
    
  return(list(xlim = xlim, ylim = ylim))
}

# Function to impute missing values. 
imputeMissingValues = function(x, impute.scale, impute.value) {
  na.index = which(is.na(x))
  if (length(na.index) > 0) {
    if (class(x) == "numeric") {
      x[na.index] = max(x, na.rm = TRUE) + impute.scale * (diff(range(x, na.rm = TRUE)))
    } 
    if (class(x) == "factor") {
      levels(x) = c(levels(x), impute.value)
      x[na.index] = impute.value
    }
  }
  return(x)
}

# subset rows and cols of the opt.path and return list with data.frames for
# x and y space and the subsets.
# returns list with dataframes for x and y space, vectors for dob, type and alpha
# character vector for x and y names
getAndSubsetPlotData = function(op, iters, subset.obs, subset.vars, subset.targets,
  marked = NULL, alpha = TRUE, impute.scale = 0.05, impute.value = "missing", ...) {
  
  # extract initial information and the data from the opt.path
  x.names = colnames(getOptPathX(op))
  y.names = op$y.names
  dim.x = length(x.names)
  dim.y = length(y.names)
  iters.max = max(getOptPathDOB(op))
  
  op.x = as.data.frame(op, include.x = TRUE, include.y = FALSE,
    include.rest = FALSE, dob = 0:max(iters), eol = c(min(iters):iters.max, NA))
  op.y = as.data.frame(op, include.x = FALSE, include.y = TRUE,
    include.rest = FALSE, dob = 0:max(iters), eol = c(min(iters):iters.max, NA))
  op.rest = as.data.frame(op, include.x = FALSE, include.y = FALSE,
    include.rest = TRUE, dob = 0:max(iters), eol = c(min(iters):iters.max, NA))
  dob = getOptPathDOB(op, dob = 0:max(iters), eol = c((max(iters) + 1):iters.max, NA))
  
  # mark best point / pareto front if marked = "best"
  if (is.character(marked)) {
    if(length(y.names) == 1) {
      marked = getOptPathBestIndex(op)
    } else {
      marked = getOptPathParetoFront(op, index = TRUE)
    }
  }
  
  # make sure that only points are marked that are alive at this iteration
  marked = marked[marked <= nrow(op.x)]
  
  # set alpha and type values
  .alpha = if(alpha && max(iters) > 0)
    normalize(dob, "range", range = c(1 / (max(iters) + 1), 1)) else rep(1, length(dob))
  .type = as.factor(ifelse(dob == 0, "init", ifelse(dob == max(iters), "prop", "seq")))
  .type = factor(.type, levels = c("init", "seq", "prop", "marked"))
  if (!is.null(marked)) {
    .type[marked] = "marked"
  }
  .alpha = pmax(0.1, .alpha)
  
  # Check and calculate the subsets
  if (missing(subset.obs))
    subset.obs = 1:nrow(op.x)
  assertIntegerish(subset.obs, lower = 1, upper = getOptPathLength(op), unique = TRUE, 
    any.missing = FALSE)
  # use only indices avaible in the current iterations
  subset.obs = subset.obs[subset.obs <= nrow(op.x)]
  
  if (missing(subset.vars))
    subset.vars = x.names
  if (is.numeric(subset.vars)) {
    assertIntegerish(subset.vars, lower = 1, upper = dim.x, unique = TRUE, any.missing = FALSE)
  }
  else 
    assertSubset(subset.vars, x.names)
  
  if (missing(subset.targets))
    subset.targets = y.names
  if (is.numeric(subset.targets)) {
    assertIntegerish(subset.targets, lower = 1, upper = getOptPathLength(op), unique = TRUE, 
      any.missing = FALSE)
  }
  else
    assertSubset(subset.targets, y.names)
  
  # impute missing values - don't impute in op.rest!
  op.x = BBmisc::dapply(op.x, fun = imputeMissingValues, impute.scale = impute.scale,
    impute.value = impute.value)
  op.y = BBmisc::dapply(op.y, fun = imputeMissingValues, impute.scale = impute.scale,
    impute.value = impute.value)
  
  # now subset everything
  op.x = op.x[subset.obs, subset.vars, drop = FALSE]
  op.y = op.y[subset.obs, subset.targets, drop = FALSE]
  op.rest = op.rest[subset.obs, -(1:2), drop = FALSE]
  dob = dob[subset.obs]
  .alpha = .alpha[subset.obs]
  .type = .type[subset.obs]
  x.names = if (is.numeric(subset.vars)) x.names[subset.vars] else subset.vars
  y.names = if (is.numeric(subset.targets)) y.names[subset.targets] else subset.targets
  
  return(
    list(
      op.x = op.x,
      op.y = op.y,
      op.rest = op.rest,
      dob = dob,
      .alpha = .alpha,
      .type = .type,
      x.names = x.names,
      y.names = y.names,
      rest.names = names(op.rest)
    )
  )
}


# Helper to get cumulated exec.time
getOptPathColAtTimes = function(op, times) {
  requirePackages("plyr")
  if (!is.data.frame(op))
    op = as.data.frame(op)
  assertNumeric(op$exec.time)
  op$exec.time.sum = cumsum(op$exec.time)
  op$finished = c(rep(FALSE, times = nrow(op) - 1), TRUE)
  plyr::adply(times, 1, getOptPathColAtTime, op = op)
}

getOptPathColAtTime = function(op.df, time) {
  col = op.df[which.last(op.df$exec.time.sum <= time), ]
  if (nrow(col) == 1) col$time = time
  col
}

