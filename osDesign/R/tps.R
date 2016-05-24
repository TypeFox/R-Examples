tps <-
function(formula = formula(data), data = sys.parent(), nn0, nn1, group, contrasts = NULL, method = "PL", cohort = TRUE, alpha = 1)
{
  call <- match.call()
  m    <- match.call(expand.dots = FALSE)
  m$method <- m$contrasts <- m$nn0 <- m$nn1 <- m$group <- m$cohort <- m$
  alpha  <- NULL
  m[[1]] <- as.name("model.frame")
  m      <- eval(m, sys.parent())
  Terms  <- attr(m, "terms")
  a      <- attributes(m)
  Y      <- model.extract(m, "response")
  X      <- model.matrix(Terms, m, contrasts)

  ## Potential Errors
  if(length(nn0) != length(nn1)) 					 stop("nn0 and nn1 should be of same length")
  if(length(nn0) != length(unique(group))) stop("Number of strata defined by group should be same as length of nn0")
  if(length(group) != nrow(X))						 stop("Group and x are not compatible")
  if((any(nn1 == 0)) || (any(nn0 == 0)))   stop("Zero cell frequency at phase I")

  ## method
  imeth <- charmatch(method, c("PL", "WL", "ML"), nomatch = 0)
  methodName <- switch(imeth + 1,
                        stop("Method doesn't exist"),
                        "PL",
                        "WL",
                        "ML")
  ## data
  if(is.matrix(Y) && (ncol(Y) > 1))
  {
    case <- as.vector(Y[, 1])
    N    <- as.vector(Y[, 1] + Y[, 2])
  }
  else
  {
    case <- as.vector(Y)
    N    <- rep(1, length(case))
  }
	
	## evaluation
	z   <- call(methodName, nn0 = nn0, nn1 = nn1, x = X, N = N, case = case, group = group, cohort = cohort, alpha = alpha)
	#out <- eval(z, local = sys.parent())
	out <- eval(z)
	names(out$coef) <- dimnames(X)[[2]]
	if(!is.null(out$cove)) dimnames(out$cove) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
	if(!is.null(out$covm)) dimnames(out$covm) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
	out$method <- method
	class(out) <- "tps"
	##
	return(out)
}
