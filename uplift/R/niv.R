######################################################################
# Adjusted Net Information Value
######################################################################

niv <- function(formula, data, subset, na.action = na.pass, 
                B = 10, direction = 1, nbins = 10, continuous = 4, plotit = TRUE, ...) {
 
  if (!inherits(formula, "formula"))
    stop("uplift: Method is only for formula objects.")
  
  if (B < 2)
    stop("uplift: Number of bootstrap samples must be greater than 1.")
  
  if (!direction %in% c(1, 2)) 
    stop("uplift: Direction must be either 1 or 2.")
   
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
  nwoe.res <- vector("list", nr_vars)
  niv.res <- matrix(nrow = nr_vars, ncol = B)
  quant <- seq(0, 1, 1 / nbins)
  
  for (j in 1:(B + 1)) {
    
    if (j > 1) {
      
      sample.ind <- sample(1:n_obs, n_obs, replace = TRUE)
      xs <- x[sample.ind, ]
      ys <- y[sample.ind]
      cts <- ct[sample.ind]
      
    } else {
      
      xs <- x
      ys <- y
      cts <- ct
      
      }
  
    for (i in 1:nr_vars) {
    
      xi <- xs[, i]
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
    
      xi.sum <- cbind(ct1.y1 = tapply(ys[cts==1], xi.c[cts==1], sum),
                      ct1.y0 = tapply(ys[cts==1]==0, xi.c[cts==1], sum),
                      ct0.y1 = tapply(ys[cts==0], xi.c[cts==0], sum),
                      ct0.y0 = tapply(ys[cts==0]==0, xi.c[cts==0], sum))
    
      ### set probability threshold to prevent division by zero y log calculation 
      crit.value <- 1E-11 
      xi.sum[is.na(xi.sum) | xi.sum == 0] <- crit.value
    
      xi.prop <- prop.table(xi.sum, 2)
    
      if (direction == 1) {
        xi.res <- cbind(xi.prop, 
                        ct1.woe = log(xi.prop[, 1] / xi.prop[, 2]), 
                        ct0.woe = log(xi.prop[, 3] / xi.prop[, 4]))
        xi.res <- cbind(xi.res, nwoe = xi.res[, 5] - xi.res[, 6])
        xi.niv <- sum((xi.res[, 1] * xi.res[, 4] - xi.res[, 2] * xi.res[, 3]) * xi.res[, 7])
      } else {
        xi.res <- cbind(xi.prop, 
                        ct1.woe = log(xi.prop[, 1] / xi.prop[, 2]), 
                        ct0.woe = log(xi.prop[, 3] / xi.prop[, 4]))
        xi.res <- cbind(xi.res, nwoe = xi.res[, 6] - xi.res[, 5])
        xi.niv <- sum((xi.res[, 2] * xi.res[, 3] - xi.res[, 1] * xi.res[, 4]) * xi.res[, 7])
      }
      
      xi.res <- round(xi.res, 4)
      xi.niv <- round(xi.niv,4)
      
      if (j == 1) nwoe.res[[i]] <- xi.res else # the nwoe is based on the entire sample
        niv.res[i, j-1] <- xi.niv   
    }
 }
  
    niv.res2 <- cbind(niv = apply(niv.res, 1, mean), 
                      penalty = apply(niv.res, 1, sd) * 1 / sqrt(B))
    niv.res2 <- cbind(niv.res2, adj_niv = niv.res2[, 1] - niv.res2[, 2])
                     
    res <- list(niv_val = 100 * round(niv.res2, 6), nwoe = nwoe.res)
    rownames(res[[1]]) <- var_names
    names(res[[2]]) <- var_names
  
    if (plotit) {
      res.plot <- res[[1]][, 1][order(res[[1]][, 1], decreasing = TRUE)]
      barplot(res.plot,
              horiz = TRUE,
              col = "blue",
              xlab = "Adjusted Net Information Value", ...)
   }
  
    return(res)
  
}
  
### END CODE