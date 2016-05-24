predict.aodml <- function(object, ..., type = c("link", "response"), se.fit = FALSE, newdata = NULL) {
  type <- match.arg(type)

## Change by RL, 28/05/2013
  
##	dat <- object$dat                     ## not needed here  
##	phi.f <- object$phi.formula           ## not needed here
##	X.b <- object$X.b     ## not needed here
##	X.phi <- object$X.phi ## not needed here

  ## old syntax: mu.f <- object$formula
  ## new syntax: get right-hand side fixed-effect formula
    mu.f <- object$formula[c(1, 3)]

  ## Some changes here too: add an "else" statement  
	if(!is.null(newdata)) {
	  X <- model.matrix(mu.f, data = newdata)
		mf <- model.frame(mu.f, data = newdata)
  	offset <- model.offset(mf)
	  }
  else {
    X <- object$X.b
    offset <- object$offset
  }
	
  b <- object$b
  
  # predict eta ("link")
	nu <- as.vector(X %*% b)
  eta <- if(is.null(offset)) nu else nu + offset
  varb <- vcov(object)
  
  vareta <- X %*% varb %*% t(X)
  se.eta <- as.vector(sqrt(diag(vareta)))
	
  # predict mu ("response")
	mu <- invlink(eta, type = object$link)
  J <- switch(object$link,
              cloglog = diag(-(1 - mu) * log(1 - mu), nrow = length(mu)),
              log = diag(mu, nrow = length(mu)),
              logit = diag(mu * (1 - mu), nrow = length(mu)),
              probit = diag(dnorm(eta), nrow = length(mu)))
  varmu <- J %*% vareta %*% J
  se.mu <- as.vector(sqrt(diag(varmu)))
  
	if(!se.fit)
    res <- switch(type,
                  response = mu,
                  link = eta)
  else
    res <- switch(type,
                  response = list(fit = mu, se.fit = se.mu),
                  link = list(fit = eta, se.fit = se.eta))
	res

}

predict.aodql <- function(object, ...) predict(object$fm, ...)
