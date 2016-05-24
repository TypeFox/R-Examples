##############################################################################
# Parse the data argument, perform basic checks, form model matrices, 
# infer parameter names and model size.
#
# Args:
#   data       - a list or a data.frame given as arg to hbglm() 
#   parsed.fm  - the parsed formula used to estimate the model 
#   parsed.fixed.fm   - the parsed fixed effects formula for the model 
#   scale.data - if TRUE then center and scale data in each col
#   predict    - TRUE if the data is for prediction - for prediction upper 
#                level data isn't needed & all groups needn't be present
#
# Returns:
#   A named list object. See the return object for details.
##############################################################################
model.hbglm <- function(data, parsed.fm, parsed.fixed.fm, 
                        scale.data = FALSE, predict = FALSE)
{
    grpID.col <- parsed.fm$grpID.col
    #grpID.col <- ifelse(predict, grpID.col, parsed.fm$grpID.col)

    # Do we have regression in the upper level? 
    need.upper.data <- parsed.fm$pooled &&
                      (parsed.fm$intercept[2] || !is.null(parsed.fm$upper.cov))
    # Do we have fixed.effects?
    have.fixed <- !is.null(parsed.fixed.fm)

    # Parse the data argument and perform basic checks
    if (is.data.frame(data)) {   # data is a single table
        num.obs <- nrow(data)
        if (is.null(data[[grpID.col]]))
            stop("Arg data (data.frame object) doesn't have 'grpID' column")
        grp.labels <- unique(data[[grpID.col]])
        num.grp <- length(grp.labels)
        if (!all(parsed.fm$var.names %in% colnames(data)))
            stop("Some random eff or upper level covariates not found in data")
        df.rand <- data[ , parsed.fm$lower.cov, drop = FALSE]
        df.rand[[parsed.fm$response]] <- data[[parsed.fm$response]]
        df.rand[[grpID.col]] <- data[[grpID.col]]
        if (have.fixed) {        # Make a table for fixed effects
          if (!all(parsed.fixed.fm$fixed.cov %in% colnames(data)))
              stop("Some fixed effects covariates not found in data")
          df.fixed <- data[ , parsed.fixed.fm$fixed.cov, drop = FALSE]
          df.fixed[[parsed.fixed.fm$response]] <- data[[parsed.fm$response]]
        } else df.fixed <- NULL
        if (need.upper.data && !predict) {
          # Form upper level data.frame & condense to remove redundant rows
          df.upper <- data[ , parsed.fm$upper.cov, drop = FALSE]
          df.upper[[grpID.col]] <- data[[grpID.col]]
          keep.rows <- sapply(grp.labels,
                              function(x) which(x == data[[grpID.col]])[1])
          df.upper <- df.upper[keep.rows, , drop=FALSE]
        }
    }
    else if (is.list(data)) {   # data is a list of 2 tables
        df.lower <- data$lower
        if (!is.data.frame(df.lower))
            stop("Arg data$lower isn't a data.frame object")
        if (is.null(df.lower) || is.null(df.lower[[grpID.col]]))
            stop("Lower level data missing OR doesn't have 'grpID' column")
        num.obs <- nrow(df.lower)
        grp.labels <- unique(df.lower[[grpID.col]])
        num.grp <- length(grp.labels)
        # Make a table for random effects
        if (!all(parsed.fm$lower.cov %in% colnames(df.lower)))
            stop("Some random eff covariates not found in data$lower")
        if (!parsed.fm$response %in% colnames(df.lower))
            stop("Reponse variable must be part of data$lower")
        df.rand <- df.lower[ , parsed.fm$lower.cov, drop = FALSE]
        df.rand[[parsed.fm$response]] <- df.lower[[parsed.fm$response]]
        df.rand[[grpID.col]] <- df.lower[[grpID.col]]
        # Make a fixed effects table
        if (have.fixed) {
          if (!all(parsed.fixed.fm$fixed.cov %in% colnames(df.lower)))
              stop("Some fixed effects covariates not found in data$lower")
          df.fixed <- df.lower[ , parsed.fixed.fm$fixed.cov, drop = FALSE]
          df.fixed[[parsed.fixed.fm$response]]<- df.lower[[parsed.fm$response]]
        } else df.fixed <- NULL
        # Upper level table
        df.upper <- data$upper
        if (need.upper.data && !predict) {
          if (!all(parsed.fm$upper.cov %in% colnames(df.upper)))
              stop("Some upper level covariates not found in data$upper")
          if (is.null(df.upper[[grpID.col]]))
              stop("Missing grpID.col in upper level data")
          # Ensure df.upper is ordered according to 'grp.labels'
          if (!all(grp.labels %in% df.upper[[grpID.col]]))
              stop("Certain groups in lower level have no data in upper level")
          ordering <- sapply(grp.labels,
                             function(x) which(x == df.upper[[grpID.col]])[1])
          df.upper <- df.upper[ordering, , drop=FALSE]
        }
    }
    else stop("Invalid arg 'data': must be a named list or a data.frame")
    ## At this point we have:
    ##    df.rand (with 'response' & grpID.col cols)
    ##    df.fixed (with 'response')
    ##    df.upper (ordered by 'grp.labels' which is unique(data[[grpID.col]]))

    # Insert dummy response column in upper level data.frame (for model.matrix)
    if (need.upper.data && !predict)
        df.upper[[parsed.fm$response]] <- rep(1, nrow(df.upper))
    else df.upper <- NULL

    # Sort df.rand so that rows of the same group are together
    ordering <- sapply(unique(df.rand[[grpID.col]]),
                       function(x) which(x == df.rand[[grpID.col]]))
    df.rand <- df.rand[as.vector(unlist(ordering)), , drop = FALSE]

    # Make model matrices
    fm.rand  <- formula(parsed.fm$formula, lhs=1, rhs=1)
    if (!parsed.fm$intercept[1]) 
        fm.rand <- formula(paste(deparse(fm.rand), "- 1"))
    mat.rand <- model.matrix(fm.rand, df.rand)

    mat.fixed <- if (have.fixed) {
        ff <- parsed.fixed.fm$formula
        if (!parsed.fixed.fm$intercept) 
            ff <- formula(paste(deparse(ff), "- 1"))
        model.matrix(ff, df.fixed)
    } else NULL

    mat.upper <- if (!is.null(df.upper)) {
        fm.upper  <- formula(parsed.fm$formula, lhs=1, rhs=2)
        model.matrix(fm.upper, df.upper)
    } else NULL

    # Split random effects data by groups into a list of lists
    # Each group has a sub-list indexed by an integer from 1 to nump.grp
    # A group's sublist has 2 labels: X and y (data matrix and response)
    data.split <- function(gid) {
        grp <- unique(df.rand[[grpID.col]])[gid]
        grp.rows <- which(grp == df.rand[[grpID.col]])
        grp.data <- list(X = mat.rand[grp.rows, , drop = FALSE], 
                         y = df.rand[[parsed.fm$response]][grp.rows])
        return(grp.data)
    }
    mat.rand.split <- lapply(1:num.grp, data.split)

    # Scaling & centering data (done to Z & group wise to X)
    if (scale.data && FALSE) {  ## TODO
        scaling.mat <- matrix(rep(1, num.grp * length(colnames(mat.rand))), 
                              nrow=num.grp)
        for (j in 1:num.grp) {
            scaled <- scale.cols(mat.rand.split[[j]]$X)
            mat.rand.split[[j]]$X <- scaled$data
            scaling.mat[j, ] <- scaled$scale.facts[, 2]/scaled$scale.facts[, 1]
            mat.rand.split[[j]]$scale.facts <- scaled$scale.facts
            mat.rand.split[[j]]$intercept <- scaled$intercept.col
        }
        attr(mat.rand.split, "scaling.mat") <- scaling.mat
        if (F && parsed.fm$pooled) {  # model has upper level
            scaledZ <- scale.cols(mat.upper)
            mat.upper <- scaledZ$data
            attr(mat.upper, "scale.facts") <- scaledZ$scale.facts  
            attr(mat.upper, "intercept")   <- scaledZ$intercept.col
        }
    }

    model <- list(
      # Parameter counts
      N = num.obs,                     # num lower level observed data
      J = num.grp,                     # num groups  
      K = length(colnames(mat.rand)),  # num rand eff covariates = num pools
      M = if(have.fixed) length(colnames(mat.fixed)) else 0,  # num fixed eff
      L = ifelse(need.upper.data, length(colnames(mat.upper)),
                 ifelse(parsed.fm$pooled, 1 , 0)), # num upper level vars/pool
      # Coefficient & group names
      rand.cov = colnames(mat.rand),   # name of random effects covariates
      fixed.cov = if(have.fixed) colnames(mat.fixed) else NULL,
      upper.cov = if(need.upper.data) colnames(mat.upper) else NULL,
      grp.labels = unique(df.lower[[grpID.col]]),  # lower/upper level ordering
      grp.indx = df.rand[[grpID.col]], # group ID of each row in mat.rand
      # Model information
      has.fixed = have.fixed,                   # true if fixed eff present
      has.upper.level = parsed.fm$pooled,       # FALSE if unpooled
      upper.reg = need.upper.data,              # TRUE if upper lev regression
      scaled.data = scale.data,           # TRUE if data is centered & scaled
      # Data matrices
      mat.rand.split = mat.rand.split,    # rand effects data as list of lists
      mat.fixed = mat.fixed,                    # fixed effects data 
      mat.upper = mat.upper,                    # upper level model matrix
      response = df.lower[[parsed.fm$response]] # lower level response vector
    )
    # Consistency check
    if (!predict && need.upper.data && model$J != nrow(df.upper))  
        stop("Mismatch between number of groups and upper level data rows")
    return(model)
}

##############################################################################
# Function to center and scale columns of a matrix. 
# To intercept cols 0 centering is applied & they are scaled to 1. 
# If no intercept column in data the no column is centered.
##############################################################################
scale.cols <- function(mat, zero.thresh = 1e-8) {
    center <- apply(mat, 2, mean)
    stddev <- if (nrow(mat) == 1) rep(1, ncol(mat)) else apply(mat, 2, sd)
center <- rep(0, ncol(mat)) 
#stddev <- if(nrow(mat) == 1) rep(1, ncol(mat)) else c(sd(mat[,1]), 1) 
stddev <- rep(1, ncol(mat)) 
    intercept.cols <- stddev <= zero.thresh
    if (!any(intercept.cols)) center <- rep(0, ncol(mat))
    else  # no centering for intercept cols
        center[intercept.cols] <- 0 
    for (j in 1:ncol(mat)) {
        if (intercept.cols[j])   # scale intercept cols to 1
            stddev[j] <- ifelse(mat[1, j] > zero.thresh, mat[1, j], 1)
            mat[ , j] <- (mat[ , j] - center[j]) / stddev[j]
    }
    scale.facts <- matrix(cbind(stddev, center), ncol = 2)
    rownames(scale.facts) <- colnames(mat)
    colnames(scale.facts) <- c("scale", "center")
    intercept.index <- NA
    if (any(intercept.cols)) intercept.index <- which(intercept.cols)[1]
    return(list(data = mat, 
                scale.facts = scale.facts,
                intercept.col = intercept.index))
}

##############################################################################
# Function to get internal variable names
#### ##########################################################################
hbglm.varnames <- function(formula, data, grpID.col = "grpID")
{
    df <- data[1, , drop=FALSE]
    model <- model.hbglm(df, parseFormula(formula), grpID.col)

    varnames <- list(
        betacol  = model$lower.cov,
        betarow  = model$grp.labels,
        thetacol = model$lower.cov,
        thetarow = if (model$L > 1) model$upper.cov else "theta", 
        tau      = model$grp.labels,
        sigma    = model$lower.cov
    )
    return(varnames)
}
