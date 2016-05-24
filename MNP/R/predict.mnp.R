predict.mnp <- function(object, newdata = NULL, newdraw = NULL, n.draws = 1,
                        type = c("prob", "choice", "order"),
                        verbose = FALSE, ...){

  if (NA %in% match(type, c("prob", "choice", "order")))
    stop("Invalid input for `type'.")
  if (n.draws < 1)
    stop("Invalid input for `n.draws'.")

  p <- object$n.alt
  if (is.null(newdraw)) 
    param <- object$param
  else
    param <- newdraw
  n.cov <- ncol(coef(object))
  coef <- param[,1:n.cov]
  n.mcmc <- nrow(coef)
  cov <- param[,(n.cov+1):ncol(param)]
  
  ## get X matrix ready
  if (is.null(newdata)) 
    x <- object$x
  else {
    call <- object$call
    x <- xmatrix.mnp(as.formula(call$formula), data = newdata,
                     choiceX = call$choiceX,
                     cXnames = eval(call$cXnames),
                     base = object$base, n.dim = p-1,
                     lev = object$alt, MoP = is.matrix(object$y),
                     verbose = FALSE, extra = FALSE)
    if (nrow(x) > 1) 
      x <- as.matrix(x[apply(is.na(x), 1, sum)==0,])
    else if (sum(is.na(x))>0)
      stop("Invalid input for `newdata'.")
  }
  n.obs <- nrow(x)
  if (verbose) {
    if (n.obs == 1)
      cat("There is one observation to predict. Please wait...\n")
    else
      cat("There are", n.obs, "observations to predict. Please wait...\n")
  }

  alt <- object$alt
  if (object$base != alt[1]) 
    alt <- c(object$base, alt[1:(length(alt)-1)])

  res <- .C("predict", as.double(x), as.integer(n.obs),
            as.double(coef), as.double(cov), as.integer(p-1),
            as.integer(n.cov), as.integer(n.mcmc),
            as.integer(n.draws), as.integer(verbose),
            prob = if (n.draws > 1) double(n.obs*n.mcmc*p)
                   else double(n.obs*p),
            choice = double(n.obs*n.mcmc),
            order = double(n.obs*n.mcmc*p),
            PACKAGE = "MNP")

  ans <- list()
  if ("choice" %in% type)
    ans$y <- matrix(res$choice, ncol=n.mcmc, nrow = n.obs,
                    byrow = TRUE, dimnames = list(rownames(x), 1:n.mcmc))
  else
    ans$y <- NULL
  if ("order" %in% type)
    ans$o <- aperm(array(res$order, dim=c(p, n.mcmc, n.obs),
                         dimnames = list(alt, 1:n.mcmc,
                           rownames(x))), c(3,1,2))
  else
    ans$o <- NULL
  if ("prob" %in% type)
    if (n.draws > 1)
      ans$p <- aperm(array(res$prob, dim = c(p, n.mcmc, n.obs),
                           dimnames = list(alt, 1:n.mcmc,
                             rownames(x))), c(3,1,2))
    else
      ans$p <- matrix(res$prob, nrow = n.obs, ncol = p, byrow = TRUE,
                      dimnames = list(rownames(x), alt))
  else
    ans$p <- NULL

  ans$x <- x
  class(ans) <- "predict.mnp"
  return(ans)
}
