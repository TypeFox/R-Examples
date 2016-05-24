######################################################################
# Response Variable Transformation for Uplift (RVTU)
######################################################################

rvtu <- function(formula, data, subset, na.action = na.pass, 
                    method = c("undersample", "oversample", "weights", "none")) {
 
  if (!inherits(formula, "formula"))
    stop("uplift: Method is only for formula objects.")
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("formula", "data", "subset", "na.action"), 
                names(mf), 0)
  method.exist <- match(c("method"), 
                        names(mf), 0)
  if (method.exist == 0)
    stop("uplift: missing method in the function call")
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
  
  ### check sampling method
  samp.method <- charmatch(tolower(method), 
                          c("undersample", "oversample", "weights", "none"))
  if (is.na(samp.method))
    stop("uplift: method must be one of 'undersample', 'oversample', 'weights' or 'none'. Aborting...")
   
  ### Proportion and counts treated
  ct.tab <- table(ct)
  n.ct0 <- ct.tab[1]
  n.ct1 <- ct.tab[2]
  prop.ct1 <- prop.table(ct.tab)[2]
  
  #if (n.ct0 == n.ct1 & method != "none")
  #  stop("uplift: method must be 'none' when number of control and treatment observations are equal. Aborting...")
  
  #if (n.ct0 != n.ct1 & method == "none")
  #  stop("uplift: method must be other than 'none' when number of control and treatment observations are not equal. Aborting...")
  
  ### Define transformed response variable
  z <- ifelse((ct == 1 & y == 1) | (ct == 0 & y == 0), 1, 0)
  dframe <- cbind(x, ct, y, z)
  dframe.s <- split(dframe, as.factor(ct))
  
  
  if (method == "none") dframe.out <- dframe
  
  if (method == "undersample")  {   
    if (n.ct1 > n.ct0) {
      ind <- sample(1:n.ct1, n.ct0, replace = FALSE)
      dframe.ct1 <- dframe.s[[2]][ind, ]
      dframe.out <- rbind(dframe.s[[1]], dframe.ct1)
    } else {
      ind <- sample(1:n.ct0, n.ct1, replace = FALSE)
      dframe.ct0 <- dframe.s[[1]][ind, ]
      dframe.out <- rbind(dframe.ct0, dframe.s[[2]])     
    }
  }  
      
  if (method == "oversample")  {      
    if (n.ct1 > n.ct0) {
      ind <- sample(1:n.ct0, n.ct1, replace = TRUE)
      dframe.ct0 <- dframe.s[[1]][ind, ]
      dframe.out <- rbind(dframe.ct0, dframe.s[[2]])
    } else {
      ind <- sample(1:n.ct1, n.ct0, replace = TRUE)
      dframe.ct1 <- dframe.s[[2]][ind, ]
      dframe.out <- rbind(dframe.s[[1]], dframe.ct1)       
    }
  }  
  
  if (method == "weights") {
     dframe$w = ifelse(dframe$ct == 1, 1 - prop.ct1, prop.ct1)
     dframe.out <- dframe
  }                  
  
  return(dframe.out)
}
  
### END CODE