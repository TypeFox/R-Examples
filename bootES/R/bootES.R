## Daniel Gerlanc and Kris Kirby (2010-2012)
## High-level function for bootstrap analyses using the 'boot' package

bootES <- function(data, R=2000, data.col=NULL, group.col=NULL, block.col=NULL,
                   effect.type=c("unstandardized", "cohens.d", "hedges.g",
                     "cohens.d.sigma", "r", "akp.robust.d"),
                   contrast=NULL,
                   slope.levels=NULL,
                   glass.control=NULL,
                   scale.weights=TRUE,
                   ci.type=c("bca", "norm", "basic", "perc", "stud", "none"),
                   ci.conf=0.95,
                   plot=FALSE,
                   ...) {

  
  ## Performs different variants of bootstrap analyses for calculating
  ## effect sizes.
  ##
  ## Args:
  ##   data          : a vector or data frame containing the one or more
  ##                   columns of values (required), group labels (optional),
  ##                   and contrast (optional) for each sample
  ##   R             : the number of bootstrap 'repetitions' to perform
  ##   data.col      : The column in 'data' containing the sample values
  ##   group.col     : The column in 'data' containing the grouping info
  ##   effect.type   : The type of standardization to perform
  ##   contrast      : A named vector specifying the lambdas for different
  ##                   groups in 'data'
  ##   slope.levels  : A named vector specifying the levels for different
  ##                   groups in 'data'
  ##   glass.control : The group for which the standard deviation should
  ##                   be used, eg. "glass.control='A'"
  ##   scale.weights : TRUE/FALSE, scale the lambdas to [-1, 1]
  ##   ci.type       : The type of confidence interval to generate
  ##                   (see 'boot.ci')
  ##   ci.conf       : The confidence level of the interval
  ##   plot          : Generate the plot?
  ##   ...           : additional arguments passed to 'boot.ci'
  ##
  ## Returns:
  ##   An object of class 'bootES' and 'boot'
  ##
  ## Details: If 'R' is not a whole number, it will be round down to the nearest
  ## whole number.  * stat=cor: 'data' must be a two-column data frame, where
  ## each of the columns is numeric

  ## Error handling
  effect.type = match.arg(effect.type)
  ci.type     = match.arg(ci.type)
  
  ## Checks on 'data'.
  if (!(is.data.frame(data) || is.numeric(data)))
    stop("'data' must be a data.frame or numeric vector.")

  ## If 'data' has been passed as a numeric vector, save it as a data.frame
  ## with a data.col='scores'
  if (is.numeric(data)) {
      data      = data.frame(scores=data, row.names=NULL)
      data.col = "scores"
  }
  
  if (!nrow(data) > 0)
    stop("'data' contains no records!")
  
  ## Checks on 'R'.
  R = as.integer(R)
  r.is.valid = (length(R) == 1) && is.numeric(R) && R > 0
  if (!r.is.valid)
    stop("R must be an integer of length 1 and greater than 0")
  
  ## Check and extract 'data.col'.
  if (!is.null(data.col)) {
    if (!is.character(data.col))
      stop("'data.col' must be a character vector.")
  
    if (!data.col %in% colnames(data))
      stop("'data.col' missing from 'data'")

    vals = data[[data.col]]
  }
  
  ## Check and extract 'group.col'.
  grps = NULL
  if (!is.null(group.col)) {
    if (!is.character(group.col)) 
      stop("'group.col' must be a character vector.")
    
    if (!group.col %in% colnames(data))
      stop("'group.col' missing from 'data'")
    
    grps = as.factor(data[[group.col]])
  }  

  ## Check and extract 'block.col'.
  blocks = NULL
  if (!is.null(block.col)) {
    if (!is.character(block.col)) 
      stop("'block.col' must be a character vector.")
    
    if (!block.col %in% colnames(data))
      stop("'block.col' missing from 'data'")
    
    blocks = as.factor(data[[block.col]])
  }
  
  ## Checks on scale.weights
  if (!is.logical(scale.weights) || length(scale.weights) != 1)
    stop("'scale.weights' must be a logical vector of length 1.")

  ## Assert that ... arguments are valid 'bootES' arguments
  dot.args = list(...)
  if (length(dot.args)) {
    dot.names = names(dot.args)
    valid.args = unique(names(c(formals(boot), formals(boot.ci))))
    invalid.args = setdiff(dot.names, valid.args)
    if (length(invalid.args)) {
      msg = sprintf("Invalid arguments to bootES: %s",
                    paste0(sQuote(invalid.args), collapse=", "))
      stop(msg)
    }
    
    if ('type' %in% dot.names) {
      stop("'type' must be specified as 'ci.type'")
    }
    
  }
  
  ## Process arguments for slope calculations
  lmbds = NULL
  lmbds.orig = contrast
  slope.levels.orig = slope.levels
  if (!is.null(slope.levels))
    effect.type = "slope"
  if (identical(effect.type, "slope")) {
    if (is.null(vals))
      stop("Invalid 'data.col'.")

    ## Process the column containing slope levels
    if (is.character(slope.levels)) {
      if (!slope.levels %in% names(data))
        stop("The column '", slope.levels, "' named in 'slope.levels'",
             "does not exist in 'data'.")

      group.col = slope.levels
      grps = as.factor(data[[slope.levels]])
      unique.levels = sort(unique(data[[slope.levels]]))
      slope.levels  = structure(unique.levels, names=unique.levels)
    }
    
    if (is.null(grps))
      stop("Invalid 'group.col'")

    if (!is.null(contrast))
      stop("Cannot specify 'contrast' and 'slope.levels'")
    
    invalid.levels <- !is.numeric(slope.levels) || is.null(names(slope.levels))
    if (invalid.levels)
      stop("'slope.levels' must be a named numeric vector.")
    
    lmbds <- calcSlopeLambdas(slope.levels)

    ## Assert that there are no NA slopes, then subset 'data' to the groups for
    ## which slope.levels were provided.
    lmbds.exist   = names(lmbds) %in% grps
    missing.lmbds = names(lmbds)[!lmbds.exist]
        
    if (length(missing.lmbds))
      stop(paste("'", missing.lmbds, "'", sep="", collapse=", "),
           " is/are not valid groups.")
    
    boot.groups = unique(names(lmbds))
    data  = data[data[[group.col]] %in% boot.groups, ]

    grps = as.factor(data[[group.col]])
    vals = data[[data.col]]
  }

  
  ## Check and extract contrast.  
  if (!is.null(contrast)) {            
    use.default.contrast = is.character(contrast) && length(contrast) == 2
    if (use.default.contrast)
      contrast = structure(c(-1, 1), names=contrast)

    if (is.null(names(contrast)))
      stop("'contrast' must be a named vector")

    if (is.null(group.col))
      stop("Must specify a 'group.col' when providing a 'contrast' argument.")
    
    lmbds = contrast

    ## Assert that contrast sum to 0
    if (!isTRUE(all.equal(sum(lmbds), 0, tol=1e-2)))
      stop("'contrast' must sum to 0.")
    
    ## Scale contrast if specified and not using the slope function
    scale.lambdas = effect.type != 'unstandardized' || scale.weights
    if (scale.lambdas)
      lmbds = scaleLambdasBySide(lmbds)    
    
    ## Assert that there are no NA lambdas, then subset 'data' to the groups for
    ## which contrast were provided.
    lmbds.exist   = names(lmbds) %in% grps
    missing.lmbds = names(lmbds)[!lmbds.exist]
        
    if (length(missing.lmbds))
      stop(paste("'", missing.lmbds, "'", sep="", collapse=", "),
           " is/are not valid groups.")
    
    boot.groups = unique(names(lmbds))
    data  = data[data[[group.col]] %in% boot.groups, ]

    vals = data[[data.col]]
    grps = as.factor(data[[group.col]])
  }

  ## Determine the 'stat' based on the passed arguments
  stat = determineStat(data, data.col, grps, effect.type, lmbds)
  
  ## Error handling for the stat='cor'
  if (stat == "cor") {
    is.valid = length(data) == 2 && all(sapply(data, is.numeric))
    if (!is.valid)
      stop("'data' must be a data frame with two numeric columns.")
  }

  ## Error handling for the stat='cor.diff'
  if (stat == "cor.diff") {
    if (is.null(group.col))
      stop("You must specify a grouping column.")

    if (length(unique(grps)) != 2)
      stop("There must be only 2 distinct groups!")

    ## Assert that there are 2 numeric columns and reorder these
    ## to be the first 2 columns in the data frame
    num.col.idx = which(sapply(data, is.numeric))
    num.col.idx = num.col.idx[!names(num.col.idx) %in% group.col]    
    is.valid = length(num.col.idx) == 2
    if (!is.valid)
      stop("'data' must contain 2 numeric columns and a grouping column.")
    
    group.col.idx = match(group.col, names(data))
    data = data[, c(num.col.idx, group.col.idx)]
  }    
  
  ## Modify the groups and lambdas when a blocking column has been passed:
  grps.char = if (!is.null(grps)) as.character(grps) else NULL
  blocks.char = if (!is.null(blocks)) as.character(blocks) else NULL
  if (!is.null(grps.char) && !is.null(blocks.char)) {
    if (!is.null(glass.control))
      stop("Cannot use 'block.col' with 'glass.control'.")
    
    blocks.char = paste(grps.char, blocks.char, sep="-")
    blocks = as.factor(blocks.char)
    block.grps = split(grps.char, blocks.char, drop=TRUE)
    ## Assert that no blocks contain values from more than one group
    n.unique.grps = sapply(block.grps, function(x) length(unique(x)))
    if (isTRUE(any(n.unique.grps > 1)))
      stop('Blocks cannot contain multiple groups.')
    block.grps = sapply(block.grps, "[", 1)
    lmbd.adj = c(table(block.grps))

    ## Adjust lambdas according to the number of blocks per group
    lmbds = lmbds / lmbd.adj[names(lmbds)]
    lmbds = lmbds[block.grps]
    names(lmbds) = names(block.grps)
    grps = blocks
    grps.char = blocks.char
    ## grp.idx = split(seq_along(vals), grps.char, drop=TRUE)
    ## block.idx = split(seq_along(vals), blocks.char, drop=TRUE)
  } else if (is.null(grps.char) && !is.null(blocks.char)) {
    grps <- as.factor(blocks.char)
    grps.char <- blocks.char
    n.by.grp <- c(table(grps.char))
    if (is.null(lmbds)) {
      lmbds <- structure(rep(1 / length(n.by.grp), length(n.by.grp)),
                         names=names(n.by.grp))
    }
  }
  
  ## Error handling for 'glass.control'
  if (!is.null(glass.control)) {
    if (!is.character(glass.control) || length(glass.control) != 1)
      stop("glass.control must be a character vector of length 1.")

    if (!glass.control %in% grps)
      stop("'glass.control' is not a valid group.")
  }  

  ## Simplest Case: No groups, so we can calculate all of the stats for a single
  ## group
  n.grps = length(unique(grps))
  single.group = (is.null(grps) || n.grps == 1) && stat != "cor"
  if (single.group) {
    boot.fun <- switch(effect.type,
                       unstandardized = meanBoot,
                       r              = rMeanBoot,
                       hedges.g       = hMeanBoot,
                       cohens.d       = dMeanBoot,
                       akp.robust.d   = akpRobustD,
                       cohens.d.sigma = dSigmaMeanBoot)    
    res = boot(vals, statistic=boot.fun, R=R) 
  } else { # Two or more groups
    use.cohens.d = effect.type %in% 
      c("cohens.d", "hedges.g", "cohens.d.sigma", "akp.robust.d")
    if (use.cohens.d) {
      res = boot(vals, calcCohensD, R, stype="f", strata=grps, grps=grps,
                 contrast=lmbds, 
                 cohens.d.sigma=(effect.type == "cohens.d.sigma"),
                 glass.control=glass.control, 
                 hedges.g=(effect.type == "hedges.g"),
                 akp.robust.d=(effect.type == "akp.robust.d")
                 )
    } else if (stat == "cor") {
      res = boot(data, statistic=corBoot, R=R) 
    } else if (stat == "cor.diff") {
      res = boot(data, calcBootCorDiff, R=R, stype="f", strata=grps, grps=grps)
    } else if (stat == "slope") {      
      res = boot(vals, calcSlope, R=R, stype="f", strata=grps,
        grps=grps, lambdas=lmbds)
    } else if (stat == "mean" && effect.type == "unstandardized") {
      res = boot(vals, meanUnweightedBoot, R, stype="f", 
                 strata=grps, grps=grps)
    } else if (stat == "contrast" && effect.type == "unstandardized") {
      res = boot(vals, calcUnstandardizedMean, R, stype="f", strata=grps,
        grps=grps, lambdas=lmbds)
    } else if (stat == "contrast" || effect.type == "r") {
      res = boot(vals, calcPearsonsR, R, stype="f", strata=grps, grps=grps,
        lambdas=lmbds)
    } else {
      fmt = paste("effect.type: %s for a 'data' of class '%s' of length %d",
        "with data.col of class %s and group.col of class %s not implemented.",
        collapse="")
      msg = sprintf(msg, effect.type, class(data), length(data), class(grps),
        class(vals))
      stop(msg)
    }
  }

  ## Calculate the confidence interval
  if (ci.type != "none") {
    ci = boot.ci(res, conf=ci.conf, type=ci.type, ...)
    ci = ci[[ci.type, exact=FALSE]]

    bounds = switch(ci.type,
      norm  = ci[1, 2:3, drop=TRUE],
      ci[1, 4:5, drop=TRUE])
  } else {
    bounds = c(NA_real_, NA_real_)
  }

  if (plot)
    plot(res)

  res[["bounds"]]  = bounds
  res[["ci.type"]] = ci.type
  res[["ci.conf"]] = ci.conf
  res[["contrast"]] = lmbds.orig
  res[["contrast.scaled"]] = lmbds
  class(res) = c("bootES", "boot")
  
  return(res)
}

print.bootES <- function(x, ...) {
  ## Prints the bootstrap summary statistics for a bootES object, followed by
  ## the confidence interval if the user required its calculation.
  printTerse(x)  
}

printTerse <- function(x) {
  ## Print a terse representation of the bootES results describing the
  ## confidence level, the type of confidence interval, the statistic, and the
  ## calculated CI bounds
  stopifnot(inherits(x, "bootES"))
  cat("\n")
  
  ## Print the scaled and unscaled contrast.
  if (!is.null(x$contrast)) {
    lmbds.txt <- gsub(" +", "", formatC(x$contrast, digits=3, format="fg"))
    cat(sprintf("User-specified lambdas: (%s)\n",
                paste(lmbds.txt, collapse=", ")))
                      
  }
  if (!is.null(x$contrast.scaled)) {
    lmbds.scl.txt <- gsub(" +", "",
                      formatC(x$contrast.scaled, digits=3, format="fg"))
    cat(sprintf("Scaled lambdas: (%s)\n",
                paste(lmbds.scl.txt, collapse=", ")))
                      
  }

  ## BEGIN: Code from boot::print.boot
  index  = seq_len(ncol(x$t))
  t      = matrix(x$t[, index], nrow = nrow(x$t))
  allNA  = apply(t, 2L, function(t) all(is.na(t)))
  ind1   = index[allNA]
  index  = index[!allNA]
  t      = matrix(t[, !allNA], nrow = nrow(t))
  t0     = x$t0

  bias = apply(t, 2L, mean, na.rm=TRUE) - t0
  std.error = sqrt(apply(t, 2L, function(t.st) var(t.st[!is.na(t.st)])))
  ## END: Code from boot::print.boot
  
  nms = c("Stat", "CI (Low)", "CI (High)", "bias", "std. error")
  res = matrix(c(t0, x[["bounds"]], bias, std.error),
               nrow=1, dimnames=list(NULL, nms))
  res.orig <- res
  res = as.vector(round(res, 3))

  cat(sprintf("%.2f%% %s Confidence Interval, %d replicates\n",
              100 * x[["ci.conf"]],
              x[["ci.type"]],
              x[["R"]]))

  hdr <- data <- character()  
  cols   <- c("Stat", "CI (Low)", "CI (High)", "bias", "SE")
  for (i in seq_along(cols)) {
    col <- cols[i]
    data[i]    <- sprintf("%-6.3f", res[i])
    char.diff <- nchar(data[i]) - nchar(col)
    spaces <- paste(rep(" ", char.diff + 5L), collapse="")
    hdr[i] <- paste(col, spaces, collapse="")
    data[i] <- paste(data[i], "     ", collapse="")
  }

  hdr <- paste(hdr, collapse="")
  data <- paste(data, collapse="")
  
  ## Print results w/ 3 digits to the right of the decimal point
  cat(hdr, data, "", sep="\n")
  invisible(res.orig)
}

determineStat <- function(data,
                          data.col=NULL,
                          grps=NULL,
                          effect.type=NULL,
                          contrast=NULL) {
  ## Based on the arguments passed to bootES, determine the type of statistic to
  ## calculate
  ##
  ## Args:
  ##  data:
  ##  data.col
  ##  grps:
  ##  effect.type:
  ##  contrast:
  ## 

  do.cor = (is.data.frame(data) && ncol(data) == 2 && is.numeric(data[[1]]) && 
            is.numeric(data[[2]]))
  res = ''
  if (is.null(grps)) {
    ## Single Group
    ## - slope -> slope
    ## - data.col == NULL -> cor
    ## - r -> r
    ## - otherwise -> mean
    if (is.null(data.col) && do.cor) {
      res = 'cor'
    } else if (identical(effect.type, 'r')) {
      res = 'r'
    } else {
      res = 'mean'
    }
  } else { 
    if (identical(effect.type, "slope")) {
      res = 'slope'
    } else if (!is.null(contrast)) {
      res = 'contrast'
    } else {
      if (is.numeric(data)) {
        res = 'mean'
      } else {
        useCorDiff = is.data.frame(data) && length(unique(grps)) == 2 &&
        ncol(data) > 2
        if (useCorDiff)
          res = 'cor.diff'
        else
          stop('Within-subject difference between correlations not implemented.')
      }
    }
  }
  return(res)
}

summary.bootES <- function(object, ...) {
  ## BEGIN: Code from boot::print.boot
  x = object
  index  = seq_len(ncol(x$t))
  t      = matrix(x$t[, index], nrow = nrow(x$t))
  allNA  = apply(t, 2L, function(t) all(is.na(t)))
  ind1   = index[allNA]
  index  = index[!allNA]
  t      = matrix(t[, !allNA], nrow = nrow(t))
  t0     = x$t0
  
  bias = apply(t, 2L, mean, na.rm=TRUE) - t0
  std.error = sqrt(apply(t, 2L, function(t.st) var(t.st[!is.na(t.st)])))
  ## END: Code from boot::print.boot
  
  nms = c("stat", "ci.low", "ci.high", "bias", "std.error")
  res = matrix(c(t0, x[["bounds"]], bias, std.error),
               nrow=1, dimnames=list(NULL, nms))
  res
}
