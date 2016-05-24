######################################################################
# Exploratory Data Analysis Tool for Uplift Modeling
######################################################################

explore <- function(formula, data, subset, na.action = na.pass, 
                    nbins = 4, continuous = 4, direction = 1) {
 
  if (!inherits(formula, "formula"))
    stop("uplift: Method is only for formula objects.")
   
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("formula", "data", "subset", "na.action"), 
                names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  if (missing(na.action)) 
     mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  special <- "trt"
  mt <- if(missing(data)) terms(formula, special) else terms(formula, special, data = data)
  mf$formula <- mt  
  mf <- eval.parent(mf)
  Terms <- attr(mf, "terms")
  attr(Terms, "intercept") <- 0 
  trt.var <- attr(Terms, "specials")$trt
  if (length(trt.var) > 0)
    ct <- mf[, trt.var] else
      stop("uplift: Formula does not include a control/treatment variable using trt().")  
  y <- model.response(mf, "numeric")
  var_names <- attributes(Terms)$term.labels[-(trt.var-1)]
  treat.names <- levels(as.factor(ct))
  x <- model.frame(terms(reformulate(var_names)),
                    mf, na.action = na.pass)     
  if(!is.numeric(y))
    stop("uplift: The left-hand-side of the formula must be a numeric vector.")
  
  if(!is.numeric(ct))
    stop("uplift: The column containing the treatment/control variable must be a numeric vector.")
  
  ### make sure we have valid data
  if (any(is.na(as.vector(y))) |
    any(is.infinite(as.vector(y))) |
    any(is.na(as.vector(ct))) |
    any(is.infinite(as.vector(ct))))
    stop("uplift: Missing values are allowed in the predictors, but not in the response or the
          treatment/control variable.")
  
  ### check number of classes
  if (length(unique(y)) != 2)
    stop("uplift: Only binary response variables are supported.")
  
  ### check response coding
  if (!all(unique(y) %in% c(0,1)))
    stop("uplift: The response variable must be coded as 0/1.")
  
  ### check treatment coding
  if (!all(unique(ct) %in% c(0,1)))
    stop("uplift: The column containing the treatment/control variable must be coded as 0/1.")
  
  ### check number of treatments
  if (length(unique(ct)) != 2)
    stop("uplift: uplift supports only 2 treatments at the moment.")
  
  ### check length of training data 
  if (length(y) != length(ct))
    stop("uplift: Length of response and treatment/control variable must be similar.")
  
  ### check number of factor levels is no more than 32
  xlevels <- lapply(x, mylevels)
  ncat <- sapply(xlevels, length)
  maxcat <- max(ncat)
  if (maxcat > 32)
    stop("uplift: Cannot handle categorical predictors with more than 32 categories.")
 
  nr_vars <- ncol(x)
  n_obs <- length(y)
  quant <- seq(0, 1, 1 / nbins)
  eda.res <- vector("list", nr_vars)

  for (i in 1:nr_vars) {
    
    xi <- x[, i]
    miss <- is.na(xi)
    if (is.numeric(xi)) {
      xi.u <- unique(xi[!miss])
      lu <- length(xi.u)
      if (lu > continuous)      
        xi.c <- cut(xi, breaks = unique(quantile(xi, quant, names = FALSE, na.rm = TRUE)), 
                    include.lowest = TRUE)
      else xi.c <- factor(xi)
      
    } else { 
      xi.c <- xi
      }
     
    # Create separate category for missing values
    if (any(miss)) {
      xi.c <- factor(xi.c, exclude = NULL)
      levels(xi.c)[length(levels(xi.c))] <- "Missing"
    }
    
    agg1.xi <- tapply(y, list(xi.c, ct), length)
    agg1.xi[is.na(agg1.xi)] <- 0 # replace NA counts with zero
    agg2.xi <- tapply(y, list(xi.c, ct), mean)
    aggAll.xi <- cbind(agg1.xi, agg2.xi)
    colnames(aggAll.xi) <- c(paste('N(Treat=', treat.names,")", sep=''),
                             paste('Mean Resp.(Treat=', treat.names,")", sep=''))
    if (direction == 1) {
      aggAll.xi <- cbind(aggAll.xi, "Uplift" = aggAll.xi[, 4] - aggAll.xi[, 3]) } else {
        aggAll.xi <- cbind(aggAll.xi, "Uplift" = aggAll.xi[, 3] - aggAll.xi[, 4]) }
                                             
    eda.res[[i]] <- round(aggAll.xi, 4)
        
  }
  names(eda.res) <- var_names
  return(eda.res)   
}
  
### END CODE