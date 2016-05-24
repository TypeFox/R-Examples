######################################################################
# Transform predictors for detecting interactions between treatment and
# covariates (based on Tial et al.)
# "A simple method for detecting interactions between treatment and
#  a large number of covariates" (2012), Tial et al. 
######################################################################

tian_transf <- function(formula, data, subset, na.action = na.pass,
                        method = c("undersample", "oversample", "none"), 
                        standardize = TRUE, cts = FALSE)  {
 
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
  oldcontrasts <- unlist(options("contrasts"))
  if (cts)
    options(contrasts = c(unordered = "contr.none", ordered = "contr.diff"))
  x <- model.matrix(terms(reformulate(var_names)),
                    mf, contrasts) 
  options(contrasts = oldcontrasts)
  intercept <- which(colnames(x) == "(Intercept)")
  if (length(intercept > 0)) x <- x[, -intercept]
  
  ### Now do all the checks...
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
                          c("undersample", "oversample", "none"))
  if (is.na(samp.method))
    stop("uplift: method must be one of 'undersample', 'oversample' or 'none'. Aborting...")
  
  ### Proportion and counts treated
  ct.tab <- table(ct)
  n.ct0 <- ct.tab[1]
  n.ct1 <- ct.tab[2]
  prop.ct1 <- prop.table(ct.tab)[2]
  
 # if (n.ct0 == n.ct1 & method != "none")
 #  stop("uplift: method must be 'none' when number of control and treatment observations are equal. Aborting...")
  
 # if (n.ct0 != n.ct1 & method == "none")
 #  stop("uplift: method must be other than 'none' when number of control and treatment observations are not equal. Aborting...")
      
  ### Standardize predictors (followed 'lars' package)
  if (standardize) {
    n <- nrow(x)
    one <- rep(1, n)
    normx <- sqrt(drop(one %*% (x^2)))
    x <- scale(x, FALSE, normx)  # scales x
  }
  
  ### Transform predictors based on Tian en al. (2012)
  trt <- ifelse(ct == 1, 1, -1)
  w <- as.data.frame(apply(x, 2, function(y) y * trt / 2))
  
  colnames(w) <- paste("T_", colnames(x), sep = '')
  dframe <- data.frame(w, ct, y)
  dframe.s <- split(dframe, as.factor(ct))
  
  ### Apply sampling method
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
   
    return(dframe.out)

}


### END CODE


#######################################
# An alternative contrasts function for unordered factors
# Ensures symmetric treatment of all levels of a factor
# (penalized package)
#######################################
#contr.none <- function(n, contrasts) {
#  if (length(n) == 1)
#    contr.treatment(n, contrasts = n<=2)
#  else
#    contr.treatment(n, contrasts = length(unique(n))<=2)
#}

#######################################
# An alternative contrasts function for ordered factors
# Ensures use of a difference penalty for such factors
# (penalized package)
#######################################
#contr.diff <- function (n, contrasts = TRUE)
#{
#  if (is.numeric(n) && length(n) == 1) {
#    if (n > 1)
#      levs <- 1:n
#    else stop("not enough degrees of freedom to define contrasts")
#  }
#  else {
#    levs <- n
#    n <- length(n)
#  }
#  contr <- array(0, c(n, n), list(levs, paste(">=", levs, sep="")))
#  contr[outer(1:n,1:n, ">=")] <- 1
#  if (n < 2)
#    stop(gettextf("contrasts not defined for %d degrees of freedom",
#                  n - 1), domain = NA)
#  if (contrasts)
#    contr <- contr[, -1, drop = FALSE]
#  contr
#}

