# adaptive sparse linear model using Figueiredo EM algorithm for adaptive sparsity (Jeffreys prior)
# Followed example from Leisch - CreatingPackages
# TODO: According to Leisch, this example is minimalistic, production code should handle more cases such as:
#       missing values, weights, offsets, subsetting, etc.  (Ie see beginning of lm() in R)
# TODO: pick better names!
# TODO: implement plot, deviance, anova, step, vcov, extractAIC, logLik
# TODO: how many of these methods above can just call lm's implementation?

aslm <- function(x, ...) UseMethod("aslm")

aslm.default <- function(x, y, ...)
  {
    x <- as.matrix(x)
    y <- as.numeric(y)

    est <- figEM(x, y, ...)

    est$fitted.values <- as.vector(x %*% est$coefficients)
    est$residuals <- y - est$fitted.values
    # TODO better rank calculation!!!!
    est$rank <- 0 # set rank to 0 for now to avoid use of qr.lm in places such as plot.lm
    est$call <- match.call()
    # TODO how to handle formula and data etc. to support summary.aslm?
    est$data <- c(x, y)
    est$contrasts <- attr(x, "contrasts")
    class(est) <- c("aslm","lm")
    est
  }

getSparseModel <- function(object)
{
  zeros = which(object$coefficients==0)
  zeronames = names(object$coefficients[zeros])
  allnames = attr(object$terms,"term.labels")
  #zeronames = attr(object$terms,"term.labels")[zeros]
  sparse.model = ". ~ ."
  contrastnames = names(object$contrasts)
  namelens = sapply(contrastnames, nchar)
  if (length(namelens) > 0) {
    sortednames = sort(namelens, decreasing=TRUE, index.return=TRUE)
  } else {
    sortednames = c()
  }
  for (nn in zeronames) {
    if (nn %in% allnames) {
      # not a factor
      sparse.model = paste(sparse.model, " -", nn, sep="")
    } else {
      if (nn == "(Intercept)") {
        sparse.model = paste(sparse.model, " -1", sep="")
      } else if (is.null(sortednames)) {
        warning(paste(nn, " had a zero coefficient, but is not in term.labels and is not a factor.  Cannot remove it from the model."))
      } else {
        factor.found = FALSE
        for (cc in contrastnames[sortednames$ix]) {
          # go through contrastnames, starting with the longest name first
          # for the case when one variable's name is the same as another variable's name with
          # extra characters added on.
          if (substr(nn,1,nchar(cc))==cc) {
            # now make sure ALL levels with this factor are 0
            # tricky for now just handle 2 level case (since we don't have to check any other coefficients)
            # TODO handle more than two levels
            sparse.model = paste(sparse.model, " -", cc, sep="")
          }
        }
      }
    }
  }
  sparse.formula = deparse(update(formula(object), sparse.model), width.cutoff=500)
  m <- lm(sparse.formula, data=object$data)
  return(m)
}

print.aslm <- function(x, ...)
  {
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients)
  }

summary.aslm <- function(object, ...)
  {
    # TODO do we have a vcov? Yes, but it's wrong (just the identity matrix)
    # so instead of calculating things ourselves, lets feed the new sparse model
    # in to lm to get lm's summary calculations
    m <-getSparseModel(object)
    res <- summary(m)
    res$lm <- m
    res$call <- object$call
    res$formula <- formula(object)
    class(res) <- c("summary.aslm","summary.lm")
    res
  }

print.summary.aslm <- function(x, ...)
  {
    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("Original formula:\n")
    print(x$formula)
    cat("\n")
    
    cat("Sparse formula:\n")
    print(formula(x$lm))
    cat("\n")
    
    cat("Residuals:\n")
    print(summary(x$residuals))
    cat("\n")

    printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)

  }


aslm.formula <- function(formula, data=list(), na.action=na.omit, ...)
  {
    mf <- model.frame(formula=formula, data=data, na.action = na.action)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    y <- model.response(mf)

    est <- aslm.default(x, y, ...)
    est$call <- match.call()
    est$formula <- formula
    # TODO return all terms or only nonzero terms
    est$terms <- attr(mf, "terms")
    est$data <- mf

    est
  }

predict.aslm <- function(object, newdata=NULL, ...)
  {
    if (is.null(newdata))
      y <- fitted(object)
    else {
        if (!is.null(object$formula)) {
            ## model has been fitted using formula interface
            x <- model.matrix(object$formula, newdata)
          }
        else {
          x <- newdata
        }
        y <- as.vector(x %*% coef(object))
      }
    y
  }

logLik.aslm <- function(object, ...) NextMethod()

nobs.aslm <- function(object, ...) NextMethod()

fit.ols.lm <- function(x, y)
  {
    lmf = lm.fit(x, y)
    beta = lmf$coefficients
    sigma = sqrt(sum((lmf$residuals)^2) / length(y))
    return(list(beta = beta, sigma = sigma))
  }

init.ones <- function(x,y)
  {
    n = length(y)
    beta <- as.numeric(matrix(1,dim(x)[2],1))
    sigma = sqrt(sum((y - x %*% beta)^2) / n)
    return(list(beta = beta, sigma = sigma))
  }

init.runif <- function(x,y)
  {
    n = length(y)
    beta <- as.numeric(matrix(runif(dim(x)[2]),dim(x)[2],1))
    sigma = sqrt(sum((y - x %*% beta)^2) / n)
    return(list(beta = beta, sigma = sigma))
  }

init.rnorm <- function(x,y)
  {
    n = length(y)
    beta <- as.numeric(matrix(rnorm(dim(x)[2],0,50),dim(x)[2],1))
    sigma = sqrt(sum((y - x %*% beta)^2) / n)
    return(list(beta = beta, sigma = sigma))
  }

# Figueiredo EM algorithm for adaptive sparsity (Jeffreys prior)
figEM <-
function(x, y, init = NULL, stopDiff = 1e-8, epsilon = 1e-6, a = 1)
{
  n = dim(x)[1]
  d = dim(x)[2]
  if (n < d) {
      stop(paste('number of observations or rows (',n,') is less than number of dimensions or columns (',d,')'))
  }

  if(is.null(init))
  {
    init = fit.ols.lm(x,y)
    #init = init.ones(x,y)
    #init = init.runif(x,y)
    #init = init.rnorm(x,y)
  }
  
  beta = init$beta
  # TODO ok, perhaps we can use a better estimate of variance of beta than max(abs(beta))/2
  # also, should this somehow be related to the scale of x or y?
  naidx = is.na(beta)
  betavar = max(abs(beta[!naidx])) / 2
  beta[naidx] = rnorm(length(beta[naidx]),0,betavar)
  sigma = sqrt(sum((y - x %*% beta)^2) / length(y))
  sigmaSqr = sigma^2

  diff = stopDiff + 1
  while(diff > stopDiff)
  {
    oldBeta = beta
    oldSigmaSqr = sigmaSqr
    
    # TODO, now that we are throwing out nonzero components, is it still necessary to add in epsilon for stability?
    U = diag((abs(beta) + epsilon),d,d)
    cols <- which((U>2*epsilon) & !is.na(U), arr.ind=TRUE)[,1]
    if (length(cols) > 0) {
      Unz <- U[cols,cols]
      xnz <- x[,cols]
      betanz = try(as.numeric(Unz %*% solve(diag(a*sigmaSqr,length(cols),length(cols)) + Unz %*% t(xnz) %*% xnz %*% Unz, Unz %*% t(xnz) %*% y))); if(class(betanz) == "try-error") break

      beta[cols] <- betanz
      beta[-cols] <- 0
    }
      sigmaSqr = sum((y - x %*% beta)^2 / n)

      diff = sqrt(sum((oldBeta - beta)^2) + (oldSigmaSqr - sigmaSqr)^2)
    }

  vcov <- diag(d)
  colnames(vcov) <- rownames(vcov) <- colnames(x)
  names(beta) <- colnames(x)

  # TODO return covariances?
  return(list(coefficients = beta,
              vcov = vcov,
              sigma = sqrt(sigmaSqr),
              df = n - d ))
}


