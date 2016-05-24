
stsm.sgf <- function(x, gradient = FALSE, hessian = FALSE, deriv.transPars = FALSE)
{
  call <- match.call()
  n <- length(x@diffy)
  w <- 2 * pi * seq(0, n - 1) / n

  if (x@model %in% c("local-level", "local-trend", "llm+seas", "BSM", "trend-cycle"))
    umcosw <- 1 - cos(w)

  pars <- get.pars(x, rescale = FALSE)
  allpars <- c(pars, get.cpar(x, rescale = FALSE), get.nopars(x, rescale = FALSE))

  allnames <- names(allpars)
  p <- length(pars)
  
  if (x@model %in% c("cycle", "trend-cycle"))
    varw <- switch(x@model, "cycle" = "var2", "trend-cycle" = "var4")

  part <- x@sgfc
  gr <- hes <- NULL

  switch(x@model,
    "local-level" = 
    {
      if (is.null(x@sgfc))
      {
        part <- cbind(var1 = 2 * umcosw, var2 = rep(1, n))
      } 
    },

    "local-trend" = 
    {
      if (is.null(x@sgfc))
      {
        part <- cbind(var1 = 4 * umcosw^2, var2 = 2 * umcosw, var3 = rep(1, n))
      } 
    },

    "llm+seas" = 
    {
      if (frequency(x@y) == 1)
        warning("The series 'x@y' is not a seasonal time series.")

      if (is.null(x@sgfc))
      {
        s <- frequency(x@y)
        umcosws <- 1 - cos(w*s)
        part3 <- umcosws / umcosw
        part3[1] <- s^2
        part <- cbind(var1 = 2 * umcosws, var2 = part3, var3 = 2 * umcosw)
      }
    },

    "BSM" = 
    {
      if (frequency(x@y) == 1)
        warning("The series 'x@y' is not a seasonal time series.")
      if (is.null(x@sgfc))
      {
        s <- frequency(x@y)
        umcosws <- 1 - cos(w*s)
        part3 <- umcosws / umcosw
        part3[1] <- s^2
        part <- cbind(var1 = 4 * umcosw * umcosws, var2 = 2 * umcosws, 
          var3 = part3, var4 = 4 * umcosw^2)
      }
    }
  )

  # spectral generating function

  ref <- charmatch(colnames(part), allnames)
  sgf <- drop(part %*% allpars[ref])

  # derivatives

  do.gr.hestpars <- (gradient || (!is.null(x@transPars) && hessian))
  if ((gradient || hessian) && !is.null(x@transPars))
    dtrans <- transPars(x, gradient = do.gr.hestpars, hessian = hessian)

  if (do.gr.hestpars) # hessian with transPars uses 'gr0'
  {
    pars.names <- names(pars)
    gr <- matrix(nrow = n, ncol = ncol(part))
    colnames(gr) <- colnames(part) #pars.names

    ref <- pars.names %in% colnames(part)
    if (any(ref))
    {
      id <- pars.names[ref]
      gr[,id] <- part[,id]
    }

    gr0 <- gr

    if (!is.null(x@transPars) && deriv.transPars)
      gr <- t(t(gr) * dtrans$gradient)
  }

  if (hessian)
  {

if (!x@model %in% c("cycle", "trend-cycle"))
{
    hes <- array(0, dim = c(n, ncol(part)))
    colnames(hes) <- colnames(part) #names(pars)

} else #for version in development
{
labels <- colnames(part) #names(pars)
hes <- array(0, dim = c(ncol(part), ncol(part), n))
dimnames(hes)[[1]] <- dimnames(hes)[[2]] <- labels
}

    if (!is.null(x@transPars) && deriv.transPars)
    {

if (!x@model %in% c("cycle", "trend-cycle")) {
      hes <- t(apply(gr0, MARGIN = 1, FUN = function(x, b) 
        x * b, b = diag(dtrans$hessian)))
} else
{
for (i in seq(n))
{
  tmp <- tcrossprod(dtrans$gradient)
  hes[,,i] <- hes[,,i] * tmp
  diag(hes[,,i]) <- diag(hes[,,i]) + gr0[i,] * diag(dtrans$hessian)
}

}
    }
  }

  list(sgf = sgf, gradient = gr, hessian = hes, constants = part)
}
