donner <- function(formula = NULL, response = NULL, weights = NULL, group = NULL, data, C = NULL){
  CALL <- match.call()
### syntax control
  if(missing(formula) & missing(response) | !missing(formula) & !missing(group))
    stop("You must specify either a formula or a response with a grouping variable.")
  if(!missing(formula) & !missing(response))
    stop("You cannot specify a formula AND a response.")
  if(!missing(response) & missing(group))
    stop("You must specify a grouping variable.")
### Build the response and analysed data set
## "formula" syntax
  if(!missing(formula)){
    f <- formula
    if(length(f) != 3)
      stop("The formula ", deparse(f), " is incomplete.")
    explain <- as.character(attr(terms(f), "variables"))[-(1:2)]
    if(length(explain) > 1){
      warning("The formula contains several explanatory variables (", paste(explain, collapse = ", "), ").\n",
              "Only the first one (", explain[1], ") was used.")
      }
    f <- formula(paste(f[2], "~", explain[1]))
    mf <- model.frame(f, data)
    Y <- model.response(mf)
# response is a vector of proportions: y/n ~ group, weights = n
    if(!is.matrix(Y)){
      f <- formula(paste(deparse(f), "+", deparse(substitute(weights))))
      mf <- model.frame(formula(f), data = data)
      names(mf)[ncol(mf)] <- "(weights)"
      w <- model.weights(mf)
      Y <- w * cbind(Y, 1 - Y)
      }
# otherwise: Y is a matrix: cbind(y, n - y) ~ group
    }

## "response" syntax
  else{
    f <- paste(deparse(substitute(response)), "~", deparse(substitute(group)), "+", deparse(substitute(weights)))
    f <- formula(f)
    mf <- model.frame(formula(f), data = data)
    names(mf)[ncol(mf)] <- "(weights)"
    Y <- model.response(mf)
    w <- model.weights(mf)
# response is a matrix: cbind(y, n - y) ~ group
    if(is.matrix(Y)){
      if(!missing(weights))
        warning("The response is a matrix: weights were ignored.")
      }
# response is a vector of proportions: y/n ~ group, weights = n
    else{
      if(any(abs(Y) > 1))
        stop("Proportions cannot be greater than 1.")
      Y <- w * cbind(Y, 1 - Y)
      }
    }

## at the end of this step, response is always a matrix: cbind(y, n - y)

## check response
  n <- rowSums(Y)
  if(any(Y < 0) | any(n <= 0))
    stop("Negative counts and NULL weights are not allowed.")

## Build data set
  datan <- data.frame(n = n, y = Y[, 1], groups = factor(mf[, 2]))
  
### Computations

# computation of rho and correction factors C
  groups <- sort(unique(datan$groups))
  lev <- levels(groups)
  N <- ntot <- ytot <- p <- nA <- 0
  datan$p <- rep(0, nrow(datan))
  for(i in seq(length(groups))){
    datatmp <- datan[datan$groups == lev[i], ]
    n <- datatmp$n ; y <- datatmp$y
    N[i] <- nrow(datatmp) ; ntot[i] <- sum(n) ; ytot[i] <- sum(y) ; p[i] <- ytot[i] / ntot[i]
    datan$p[datan$groups == lev[i]] <- p[i]
    nA[i] <- sum(datan$n[datan$groups == lev[i]]^2) / ntot[i]
    }
  df.SSC <- sum(N) - length(groups)
  df.SSE <- sum(ntot) - sum(N)
  MSC <- sum(datan$n * (datan$y / datan$n - datan$p)^2) / df.SSC
  MSE <- sum(datan$n * (datan$y / datan$n) * (1 - datan$y / datan$n)) / df.SSE
  K <- (sum(ntot) - sum(nA)) / (sum(N) - length(groups))
  rho <- (MSC - MSE) / (MSC + (K - 1) * MSE)
  if(is.null(C))
    C <- 1 + (nA - 1) * rho

  # results
  tab <- data.frame(groups = groups, N = N, n = ntot, y = ytot, p = p, C = C)   
  names(tab)[1] <- as.character(attr(terms(f), "variables"))[3]
  p <- sum(tab$y) / sum(tab$n)
  X2 <- sum((tab$y - tab$n * p)^2 / (tab$C * tab$n * p * (1 - p)))

  # outputs
  new(Class = "drs", CALL = CALL, tab = tab, rho = rho, X2 = X2)
  }


