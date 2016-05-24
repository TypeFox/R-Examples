
barrier.eval <- function(model, type = c("1", "2"), mu = 0.01, 
  gradient = FALSE, hessian = FALSE)
{
  #stopifnot((!is.null(model@lower) || !is.null(model@upper)))
  if (mu <= 0) {
    warning("'mu' should be a positive value or zero.",
    "'mu' was set equal to zero.")
    #barrier$mu <- 0
    return(list(barrier = 0, bl = 0, bu = 0, 
      dl1 = NULL, du1 = NULL, dl2 = NULL, du2 = NULL))
  } else # && mu > 0
  if (!is.null(model@transPars) && model@transPars == "square")
  {
    warning(paste("Parameterization", sQuote(model@transPars), 
    "does not require a barrier term."))
  }

  type <- match.arg(type)[1]
  bl <- bu <- 0
  dl1 <- du1 <- dl2 <- du2 <- NULL

  pars <- model@pars
  nmspars <- names(pars)

  lower <- model@lower
  lower <- lower[is.finite(lower)]
  upper <- model@upper
  upper <- upper[is.finite(upper)]

  if (!is.null(lower))
  {
    refl <- charmatch(names(lower), nmspars)
    bl <- switch(type, 
      "1" = -sum(log(pars[refl] - lower)), 
      "2" = sum(1 / (pars[refl] - lower)))
  }

  if (!is.null(upper))
  {
    refu <- charmatch(names(upper), nmspars)
    bu <- switch(type, 
      "1" = -sum(log(upper - pars[refu])), 
      "2" = sum(1 / (upper - pars[refu])))
  }

  if (gradient)
  {
    dl1 <- du1 <- rep(0, length(pars))
    if (!is.null(lower)) {      
      dl1[refl] <- switch(type, 
        "1" = -1 / (pars[refl] - lower), 
        "2" = -1 / (pars[refl] - lower)^2)
    }
    if (!is.null(upper)) {
      du1[refu] <- switch(type, 
        "1" = 1 / (upper - pars[refu]), 
        "2" = 1 / (upper - pars[refu])^2)
    }
  }

  if (hessian)
  {
    dl2 <- du2 <- rep(0, length(pars))
    if (!is.null(lower)) {
      dl2[refl] <- switch(type, 
        "1" = 1 / (pars[refl] - lower)^2, 
        "2" = 2 / (pars[refl] - lower)^3)
    }
    if (!is.null(upper)) {
      du2[refu] <- switch(type, 
        "1" = 1 / (upper - pars[refu])^2, 
        "2" = 2 / (upper - pars[refu])^3)
    }
  }

  list(barrier = mu * (bl + bu), bl = bl, bu = bu, 
    dl1 = dl1, du1 = du1, dl2 = dl2, du2 = du2)
}
