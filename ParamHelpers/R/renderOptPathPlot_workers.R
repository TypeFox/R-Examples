# Plot methods for renderOptPathPlot. Nearly the same interface for all plot functions
# Not all functions do have all arguments.

# @param op 
#   The optimization path
# @param .alpha [\code{numeric}]\cr
#   Vector of alpha values for the points in the plots.
# @param .type [\code{factor}]\cr
#   Vector of types of the points, factor levels are init, seq, prob and marked.
# @param dob [\code{numeric}]\cr
#   Vector of dobs
# @param log[\code{character}]\cr
#   Vector of variables to be logarithmized
# @param names [\code{character}]\cr
#   Vector of the names of the variables. Used to identify variables for the plot.
# @param short.names [\code{character}]\cr
#   Vector of the short names of the variables. This names will be printed in the plot.
#   Must be the same length as names
# @param space[\code{character}]
#   If the X-Space is plotted, space = "x", if the Y-Space is plotted, space = "y".
#   Special case 1D -> 1D also "both" is possible.
# @param iter [\code{integer(1)}]\cr
#   Current iteration.
# @param classes [\code{character(2)}]\cr
#  Classes of the variables (numeric or factor) in 2D plots.
# @param xlim, ylim [\code{numeric(2)}]\cr
#  Limits for the x and y axis respectively.
# @param colours [\code{character(4)}]\cr
#   Colours of the points/lines for the three point types init, seq, prob and marked.
# @param size [\code{numeric(1)} | NULL]\cr
#   Size of points / lines. 
# @return A ggplot object.


# Plot method for a one-dimensional numeric X- or Y-Space
# Here we use geom_density and geom_rug
# And we know both names and short.names have length 1
plot1DNum = function(op, .alpha, .type, log, names, short.names,
  space, iter, xlim, colours, ggplot.theme) {
  
  op$.alpha = .alpha
  op$.type = .type
  
  if (space == "x") {  
    title = ggplot2::ggtitle("X-Space")    
  } 
  if (space == "y") {    
    title = ggplot2::ggtitle("Y-Space")    
  }
  
  pl = ggplot2::ggplot(op, ggplot2::aes_string(x = names))
  pl = pl + ggplot2::geom_density(colour = "black")
  pl = pl + title
  pl = pl + ggplot2::xlab(short.names)
  pl = pl + ggplot2::geom_rug(ggplot2::aes_string(alpha = ".alpha", colour = ".type"), 
    sides = "b", size = 2L, data = op)
  if (names %in% log)
    pl = pl + ggplot2::coord_trans(xtrans = "log10", limx = xlim)
  else
    pl = pl + ggplot2::coord_cartesian(xlim = xlim) 
  pl = pl + ggplot2::guides(alpha = FALSE)
  pl = pl + ggplot2::scale_alpha_continuous(range = c(max(1 / (iter + 1), 0.1), 1L))
  pl = pl + ggplot2::scale_colour_manual(name = "type",
    values = c(init = colours[1L], seq = colours[2L], prop = colours[3L], marked = colours[4L]))
  pl = pl + ggplot.theme
  
  return(pl)
}


# Plot method for a one-dimensional discrete X- or Y-Space
# Here we use geom_bar
plot1DDisc = function(op, .alpha, .type, log, names, short.names,
  space, iter, ylim, colours, ggplot.theme) {
  
  op$.alpha = as.factor(.alpha)
  op$.type = .type
  
  if (space == "x") {  
    title = ggplot2::ggtitle("X-Space")    
  } 
  if (space == "y") {    
    title = ggplot2::ggtitle("Y-Space")    
  }
  
  pl = ggplot2::ggplot(op, ggplot2::aes_string(x = names[1L], fill = ".type", alpha = ".alpha"))
  pl = pl + ggplot2::geom_bar()
  pl = pl + title
  pl = pl + ggplot2::xlab(short.names)
  pl = pl + ggplot2::ylim(ylim)
  pl = pl + ggplot2::scale_alpha_discrete(range = c(max(1 / (iter + 1), 0.1), 1L))
  pl = pl + ggplot2::scale_fill_manual(name = "type",
    values = c(init = colours[1L], seq = colours[2L], prop = colours[3L], marked = colours[4L]))
  pl = pl + ggplot.theme
  pl = pl + ggplot2::guides(alpha = FALSE)
  
  return(pl)
}


# Plot method for a two-dimensional X- or Y-Space
# We use geom_point and jitter for discrete variables
# y.name: we can plot contour-lines for a singel y-variable if both x-variables
# are numeric. in this case, op.y is the data.frame containing the y.variable
plot2D = function(op, .alpha, .type, log, names, short.names, y.name = NULL, op.y = NULL,
  space, iter, classes, xlim, ylim,  colours, size, ggplot.theme) {
  
  op$.alpha = .alpha
  op$.type = .type
  
  if (space == "x") {
    title = ggplot2::ggtitle("X-Space")
  } 
  if (space == "y") {
    title = ggplot2::ggtitle("Y-Space")
  }
  if (space == "both") {
    title = ggplot2::ggtitle("X- and Y-Space")
  }
  
  factor.classes = classes == "factor"
  if (any(factor.classes)) {
    # Jitter only in the discrete directions
    pos = ggplot2::position_jitter(w = 0.1 * factor.classes[1],
      h = 0.1 * factor.classes[2])
  } else {
    pos = "identity"
  }
  
  # prepare contour plot
  if (!is.null(y.name)) {
    requirePackages(c("akima", "reshape2"), why = "renderOptPathPlot plot2D")
    fld = with(cbind(op, op.y), akima::interp(x = get(names[1L]), y = get(names[2L]), z = get(y.name)))
    df = reshape2::melt(fld$z, na.rm = TRUE)
    names(df) = c(names, y.name)
    df[[names[1L]]] = fld$x[df[[names[1L]]]]
    df[[names[2L]]] = fld$y[df[[names[2L]]]]
  }
  
  pl = ggplot2::ggplot()
  pl = pl + ggplot2::geom_point(data = op, ggplot2::aes_string( x = names[1L], y = names[2L],
    shape = ".type", colour = ".type", alpha = ".alpha"), size = size, position = pos)
  # add contour
  if (!is.null(y.name))
    pl = pl + ggplot2::stat_contour(ggplot2::aes_string(x = names[1L], y = names[2L], z = y.name), data = df)
  pl = pl + title
  pl = pl + ggplot2::xlab(short.names[1L]) + ggplot2::ylab(short.names[2L])
  pl = pl + ggplot2::guides(alpha = FALSE)
  pl = pl + ggplot2::scale_colour_manual(name = "type",
    values = c(init = colours[1L], seq = colours[2L], prop = colours[3L], marked = colours[4L]))
  pl = pl + ggplot2::scale_shape_manual(name = "type", 
    values = c(init = 15L, seq = 16L, prop = 17L, marked = 18L))
  pl = pl + ggplot2::scale_alpha_continuous(range = c(max(1 / (iter + 1), 0.1), 1L))
  pl = pl + ggplot.theme
  if (classes[1L] == "numeric") {
    if (names[1L] %in% log)
      pl = pl + ggplot2::scale_x_log10(limits = xlim)
    else
      pl = pl + ggplot2::xlim(xlim)
  }
  if (classes[2L] == "numeric") {
    if (names[2L] %in% log)
      pl = pl + ggplot2::scale_y_log10(limits = ylim)
    else
      pl = pl + ggplot2::ylim(ylim)
  }
  
  return(pl)
}

# Plot method for a multi-dimensional X- or Y-Space
# Here we make a PCP using GGally::ggparcoord
plotMultiD = function(op, .alpha, .type, log, names, short.names,
  space, iter, colours, size, scale, ggplot.theme) {
  args = list(columns = seq_along(names))
  
  # make every variable numeric and check for a log trafo
  for (var in names) {
    op[, var] = as.numeric(op[, var])
    if (var %in% log)
      op[, var] = log10(op[, var])
  }
  
  op$.alpha = .alpha
  # minimal alpha value:
  op$.type = .type
  args$data = op
  args$alphaLines = ".alpha"
  args$groupColumn = ncol(op)
  args$scale = scale
  args$mapping = ggplot2::aes_q(lwd = size)
  

  
  if (space == "x") {
    title = ggplot2::ggtitle("X-Space")
  } else {
    title = ggplot2::ggtitle("Y-Space")
  }
  pl = do.call(GGally::ggparcoord, args)
  pl = pl + ggplot2::ylab ("scaled values")
  pl = pl + ggplot2::scale_x_discrete(labels = short.names)
  pl = pl + title
  pl = pl + ggplot2::guides(alpha = FALSE, size = FALSE)
  pl = pl + ggplot2::scale_colour_manual(name = "type",
    values = c(init = colours[1L], seq = colours[2L], prop = colours[3L], marked = colours[4L]))
  pl = pl + ggplot.theme
  return(pl)
}


# Function to plot one or more numeric variables over time
# names: all corresponding variables must be numeric
# short.names: short names of the variables given by names
multiVariablesOverTime = function(op, .alpha, dob, log, names, short.names,
  space, iter, colours, ggplot.theme) {
  
  # For rest variables, we can get a NA data.frame here. In this case, no plot
  if (all(is.na(op[, names])))
    return(NULL)
  
  # allow only log trafo of all variables in this plot
  log.var = names %in% log
  if (any(log.var) && !all(log.var))
    stop("If you want to apply a log trafo in an over.time.plot, you have to apply it to every variable.")
  
  for (var in names) {
    if (!is.numeric(op[, var]))
      warning(paste("Converting variable ", var, "to numeric for over time plot."))
    op[, var] = as.numeric(op[, var])
  }
  
  op2 = op[, names]
  op2$dob = dob
  op2$.alpha = .alpha
  
  # mean over dob
  op2 = aggregate(op2, list(op2$dob), mean)[, -1]
  
  # reshape into long format
  op2 = reshape(op2, ids = row.names(op2),
    times = names, timevar = "variable",
    varying = list(names), direction = "long", v.names = c("value"))
  
  pl = ggplot2::ggplot(op2, ggplot2::aes_string(x = "dob", y = "value", group = "variable", 
    linetype = "variable"))
  pl = pl + ggplot2::geom_point()
  pl = pl + ggplot2::geom_line() 
  pl = pl + ggplot2::scale_linetype_discrete(labels = short.names)
  # For the x axis: only whole numbers as breaks
  pl = pl + ggplot2::scale_x_continuous(breaks = function(x) pretty(x, n = min(5, iter + 1)))
  
  # fixed number of decimals:
  fmt <- function(){
    function(x) format(x, nsmall = 3, scientific = FALSE)
  }
  
  if (all(log.var)) {
    pl = pl + ggplot2::scale_y_log10(labels = fmt())
  } else {
    pl = pl + ggplot2::scale_y_continuous(labels = fmt())
  }
  pl = pl + ggplot.theme
  
  return(pl)
}

# Plots One variable versus the DOB. name is the name of the variable to be plotted
oneVariableOverTime = function(op, .alpha, .type, dob, log, names, short.names, iter,
  size.points, size.lines, colours, ggplot.theme) {
  
  # For rest variables, we can get a NA data.frame here. In this case, no plot
  if (all(is.na(op[, names])))
    return(NULL)
  
  # convert factor variables to numeric
  if (!is.numeric(op[, names]))
    warning(paste("Converting variable ", names, "to numeric for over time plot."))
  op[, names] = as.numeric(op[, names])
  
  # Some data  preproc. 2 Different datasets - one for init design, one for rest
  op = cbind(op, dob = dob, .alpha = .alpha, .type = .type)
  
  init.des.inds = dob == 0
  
  op.init.des = op[init.des.inds, , drop = FALSE]
  op.seq.opt = op[!init.des.inds, , drop = FALSE]
  
  # if we want to log and all values are negative, make them positive.
  # this is special treatment for our ei.
  if (names %in% log && all(na.omit(op[, names] <= 0))) {
    op.init.des[, names] = -op.init.des[, names]
    op.seq.opt[, names] = -op.seq.opt[, names]
  }
  
  aes.points = ggplot2::aes_string(x = "dob", y = names, shape = ".type",
    colour = ".type", alpha = ".alpha")
  
  pl = ggplot2::ggplot(op, ggplot2::aes_string(x = "dob", y = names))
  # add initial design points allays with jitter in x-direction,
  # if discrete also with jitter in y-direction
  if (length(na.omit(op.init.des[, names])) > 0L) {
    if (is.numeric(op[, names]))
      pl = pl + ggplot2::geom_point(data = op.init.des, mapping = aes.points, size = size.points,
        position = ggplot2::position_jitter(height = 0.1))
    else
      pl = pl + ggplot2::geom_point(data = op.init.des, mapping = aes.points, size = size.points,
        position = ggplot2::position_jitter(height = 0.1, width = 0.1))
  }
  # add sequential points, if discrete with jitter in y-direction
  # Add jitter for discrete variable
  if (length(na.omit(op.seq.opt[, names])) > 0L) {
    if (is.numeric(op[, names]))
      pl = pl + ggplot2::geom_point(data = op.seq.opt, mapping = aes.points, size = size.points)
    else
      pl = pl + ggplot2::geom_point(data = op.seq.opt, mapping = aes.points, size = size.points,
        position = ggplot2::position_jitter(height = 0.1, width = 0.1))

    # mean data for line plot for sequential data - only for numeric vars
    # Also ylims are only useful for numeric vars
    if (is.numeric(op[, names])) {
      op.seq.means = op.seq.opt[!duplicated(op.seq.opt$dob), ]
      op.seq.means[, names] = tapply(op.seq.opt[, names], op.seq.opt[, "dob"], mean)
      pl = pl + ggplot2::geom_line(data = op.seq.means, ggplot2::aes_string(x = "dob", y = names), alpha = 0.3)
    }
  }
  
  # fixed number of decimals:
  fmt <- function(){
    function(x) format(x, nsmall = 3, scientific = FALSE)
  }
  
  if (names %in% log) {
    pl = pl + ggplot2::scale_y_log10(labels = fmt())
  } else {
    pl = pl + ggplot2::scale_y_continuous(labels = fmt())
  }
  pl = pl + ggplot2::geom_vline(xintercept = 0.5)
  pl = pl + ggplot2::guides(alpha = FALSE)
  pl = pl + ggplot2::ylab(short.names)
  pl = pl + ggplot2::scale_colour_manual(name = "type",
    values = c(init = colours[1L], seq = colours[2L], prop = colours[3L], marked = colours[4L]))
  pl = pl + ggplot2::scale_shape_manual(name = "type", 
    values = c(init = 15L, seq = 16L, prop = 17L, marked = 18L))
  
  # set range for alpha scale, so that extra variables (that may not exist in 
  # iteration 0) will have the same alpha values as all other variables.
  range = c(max(min(op$.alpha[!is.na(op[, names])]), 0.1), 1L)
  pl = pl + ggplot2::scale_alpha_continuous(range = range)
  # For the x axis: only whole numbers as breaks
  pl = pl + ggplot2::scale_x_continuous(limits = c(-0.5, NA_real_),
    breaks = function(x) pretty(x, n = min(5, iter + 1)))
  pl = pl + ggplot.theme

  return(pl)
}

