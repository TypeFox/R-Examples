create.design <- function(m, formula, nonpen  = ~ 1, data,
                          weights, subset, na.action, contrasts, env)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 30 Jun 2006, 09:46

  ## Case where some variables won't be penalized -> merge formulas,
  ## also check the environments (is this the correct way ???)

  any.nonpen <- !is.null(nonpen)
  
  ## Case where some variables won't be penalized -> merge formulas,
  ## also check the environments (is this the correct way ???)
  if(!inherits(formula, "formula") || length(formula) != 3)
    stop("Argument 'formula' is of wrong type or length")

  if(any.nonpen){
      if(!inherits(nonpen, "formula"))
        stop("Argument 'nonpen' of wrong type")

      is.gEnv <- function(e) identical(e, .GlobalEnv)

      ## Paste the two formulas together
      f <- as.formula(paste(c(deparse(formula[[2]]), "~",
                              deparse(formula[[3]]), "+",
                              deparse(nonpen[[length(nonpen)]])),
                            collapse = ""))

      ## Get the environment of the formulas
      env.formula <- environment(formula)
      env.nonpen <- environment(nonpen)

      ## If env. of 'formula' is not global, check if env. of 'nonpen' differs.
      ## If yes give warning and use env. of 'formula'
      if(!is.gEnv(env.formula)){
        environment(f) <- env.formula
        if(!is.gEnv(env.nonpen) && !identical(env.formula, env.nonpen))
          warning("'formula' and 'nonpen' have different environments. Use environment of 'formula'")
        
        environment(f) <- environment(formula)
      }else if (!is.gEnv(env.nonpen)){ ## if env. of 'nonpen' is not global
                                       ## (but env. of 'formula' is), use
                                       ## env. of 'nonpen'
        environment(f) <- env.nonpen
      }
      m$formula <- f
  }else{
    m$formula <- formula
  }
  
  m$drop.unused.levels <- TRUE
  
  ## Create model-frame
  m[[1]] <- as.name("model.frame")

  mf <- eval(m, env)

  ## Create design matrix, na.action handles the missing values, therefore
  ## weights and offset may be of shorter length (mf does this for us)
  Terms <- attr(mf, "terms") ##terms(m$formula, data = data)
  x <- model.matrix(Terms, data = mf, contrasts = contrasts)
  y <- model.response(mf)
  w <- model.weights(mf)
  off <- model.offset(mf)

  if (!is.null(w) && any(w < 0)) 
    stop("Negative weights not allowed")

  if(!is.numeric(off))
    off <- rep(0, length(y))
  if(!length(w))
    w <- rep(1, length(y))
  
  ## Handle the non-penalized coefficients
  if(any.nonpen){
    tmp <- terms(nonpen, data = data)
    co <- contrasts[attr(tmp, "term.labels")]
    
    if(length(co))
      used.co <- co[!unlist(lapply(co, is.null))]
    else
      used.co <- NULL

    ## also uses the response...to be changed
    x.nonpen <- model.matrix(tmp, data = mf, contrasts = used.co)
    matches <- match(colnames(x.nonpen), colnames(x))
  }
  index <- attr(x, "assign")
  if(any.nonpen)
    index[matches] <- NA

  list(x = x,
       y = y,
       w = w,
       off = off,
       mf = mf,
       index = index,
       Terms = Terms)
}

blockstand <- function(x, ipen.which, inotpen.which)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  4 Aug 2006, 08:50
  
  n <- nrow(x)
  x.ort <- x
  scale.pen <- list(); length(scale.pen) <- length(ipen.which)
  scale.notpen <- NULL
  
  if(length(inotpen.which) > 0){
    one <- rep(1, n)
    scale.notpen <- sqrt(drop(one %*% (x[,inotpen.which]^2)) / n)
    x.ort[,inotpen.which] <- scale(x[,inotpen.which], FALSE, scale.notpen)
  }
    
  for(j in 1:length(ipen.which)){
    ind <- ipen.which[[j]]
    decomp <- qr(x[,ind])
    if(decomp$rank < length(ind)) ## Warn if block has not full rank
      stop("Block belonging to columns ", paste(ind, collapse = ", "),
              " has not full rank! \n")
    scale.pen[[j]] <- qr.R(decomp) * 1 / sqrt(n)
    x.ort[,ind] <- qr.Q(decomp) * sqrt(n)
  }
  list(x = x.ort, scale.pen = scale.pen, scale.notpen = scale.notpen)
}

lambdamax <- function(x, ...)
  UseMethod("lambdamax")

lambdamax.formula <- function(formula, nonpen  = ~ 1, data,
                              weights, subset, na.action,
                              coef.init,
                              penscale = sqrt, model = LogReg(),
                              center = TRUE, standardize = TRUE,
                              contrasts = NULL, nlminb.opt = list(), ...)
{
  ## Purpose: Function to find the maximal value of the penalty parameter
  ##          lambda
  ## ----------------------------------------------------------------------
  ##
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 20 Apr 2006, 11:24

  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  
  ## Remove not-needed stuff to create the model-frame
  m$nonpen <- m$lambda <- m$coef.init <- m$penscale <- m$model <- m$center <-
    m$standardize <- m$contrasts <- m$... <- NULL

  l <- create.design(m, formula, nonpen, data, weights, subset, na.action,
                     contrasts, parent.frame())

  if(missing(coef.init))
    coef.init <- rep(0, ncol(l$x))
  
  lambdamax.default(l$x, y = l$y, index = l$index, weights = l$w,
                    offset = l$off, coef.init = coef.init,
                    penscale = penscale,
                    model = model, center = center, standardize = standardize,
                    nlminb.opt = nlminb.opt, ...)
}

lambdamax.default <- function(x, y, index, weights = rep(1, length(y)),
                              offset = rep(0, length(y)),
                              coef.init = rep(0, ncol(x)),
                              penscale = sqrt, model = LogReg(),
                              center = TRUE, 
                              standardize = TRUE, nlminb.opt = list(), ...)
{
  ## Purpose: Function to find the maximal value of the penalty parameter
  ##          lambda
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## X: design matrix (including intercept), already rescaled and
  ##    possibly blockwise orthonormalized.
  ## y: response vector
  ## index: vector which defines the grouping of the variables. Components
  ##        sharing the same number build a group. Non-penalized
  ##        coefficients are marked with "NA".
  ## weights: vector of observation weights.
  ## offset: vector of offset values.
  ## coef.init: initial parameter vector. Penalized groups are discarded.
  ## penscale: rescaling function to adjust the value of the penalty
  ##           parameter to the degrees of freedom of the parameter group.
  ## model: an object of class "grpl.model" implementing
  ##        the negative log-likelihood, gradient, hessian etc. See
  ##        "grpl.model" for more details.
  ## nlminb.opt: arguments to be supplied to "nlminb".
  ## ... : additional arguments to be passed to the functions defined in
  ##       model.
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 20 Apr 2006, 11:24

  any.notpen <- any(is.na(index))
  coef.npen <- coef.init[is.na(index)] ## unpenalized parameters
  
  inotpen.which <- which(is.na(index))
  
  ## Index vector of the penalized parameter groups
  ipen <- index[!is.na(index)]
  
  ## Table of degrees of freedom
  dict.pen <- sort(unique(ipen))
  ipen.tab <- table(ipen)[as.character(dict.pen)] 

  ## Indices of parameter groups
  ipen.which <- split((1:ncol(x))[!is.na(index)], ipen)

  intercept.which <- which(apply(x == 1, 2, all))
  has.intercept   <- length(intercept.which)

  if(!has.intercept & center){
    message("Couldn't find intercept. Setting center = FALSE.")
    center <- FALSE
  }

  if(length(intercept.which) > 1)
    stop("Multiple intercepts!")

  if(center){
    if(!has.intercept) ## Could be removed; already handled above
      stop("Need intercept term when using center = TRUE")

    ##if(length(intercept.which) == 1 & !is.na(index[intercept.which]))
    ##  stop("Need unpenalized intercept")

    mu.x <- apply(x[,-intercept.which], 2, mean)
    x[,-intercept.which] <- sweep(x[,-intercept.which], 2, mu.x)
  }
  
  if(standardize){
    stand        <- blockstand(x, ipen.which, inotpen.which)
    x            <- stand$x
  }

  x.npen <- x[,inotpen.which, drop = FALSE]

  helper <- function(par)
    model@nloglik(y, offset + x.npen %*% par, weights, ...)

  if(any.notpen){
    par0 <- do.call(nlminb, args = c(list(start = coef.npen,
                              objective = helper), nlminb.opt))$par
    mu0  <- model@invlink(offset + x.npen %*% par0)
  }else{
    mu0 <- model@invlink(offset)
  }
  
  ngrad0 <- model@ngradient(x, y, mu0, weights, ...)[!is.na(index)]
  
  ##gradnorms <- numeric(length(dict.pen))

  gradnorms <- c(sqrt(rowsum(ngrad0^2, group = ipen))) / penscale(ipen.tab)
                                                       
  ##for(j in 1:length(dict.pen)){
  ##  gradnorms[j] <- sqrt(crossprod(ngrad0[which(index == dict.pen[j])])) /
  ##    penscale(sum(index == dict.pen[j], na.rm = TRUE))
  ##}
  max(gradnorms)
}


