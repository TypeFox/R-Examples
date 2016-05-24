setMethod("cor.test", signature(x = "Nri"),
          function(x, y, ...)
{
  cor_apply_lh <- function(response, ...)
  {
    predictor <- get("predictor", envir= environment_apply)
    t_res <- cor.test(x = predictor, y = response, ...)
    return(c(t_res$p.value, t_res$estimate))
  }
  cor_apply_rh <- function(predictor, ...)
  {
    response <- get("response", envir= environment_apply)
    t_res <- cor.test(x = predictor, y = response, ...)
    return(c(t_res$p.value, t_res$estimate))
  }
  cor_apply_both <- function(xy, ...)
  {
    t_res <- cor.test(x = xy[1:(length(xy)/2)], y = xy[(length(xy)/2+1):length(xy)], ...)
    return(c(t_res$p.value, t_res$estimate))
  }
  
  if (class(x) == "Nri")
  {
    if (class(y) == "Nri")
    {
      nri_response <- NA
      for (i in 1:length(dim(x$nri)))
        if (dim(x$nri)[i]!=dim(y$nri)[i])
          stop("Dimensions of nri values in x and y differ")
    } else {
      nri_response <- FALSE
    }
    nri_data <- x
  } else {
    if (class(y) == "Nri")
    {
      nri_response <- TRUE
      nri_data <- y
    } else {
      stop("Could not determine which variable contains nri-values")
    }
  }
   
  if (is.na(nri_response))
  {
    xy <- new("DistMat3D", values = c(x$nri, y$nri),
              nlyr = dim(x$nri)[3] + dim(y$nri)[3],
              ncol = dim(x$nri)[1])
    res <- new("Nri", nri = xy, fwhm = x@fwhm,
                wavelength = x@wavelength,
                dimnames = list(Band_1 = x@dimnames[[1]],
                                Band_2 = x@dimnames[[1]],
                                Sample = c(x@dimnames[[1]],
                                           y@dimnames[[1]])))
    cor_data <- apply(xy, MARGIN = c(1, 2), FUN = cor_apply_both, ...)
  } else {
    if (nri_response)
    {
      res <- y
      response <- y$nri    
      environment_apply <- new.env(parent = .GlobalEnv)
      assign("predictor", x, environment_apply)

      cor_data <- apply(response, MARGIN = c(1, 2), FUN = cor_apply_lh, ...)
    } else {
      res <- x
      predictor <- x$nri
    
      environment_apply <- new.env(parent = .GlobalEnv)
      assign("response", y, environment_apply)

      cor_data <- apply(predictor, MARGIN = c(1, 2), FUN = cor_apply_rh, ...)
    }
  }
  
  final <- list(p.value = new("DistMat3D", values = cor_data[1,],
                              ncol = dim(x$nri)[1], nlyr = 1), 
                estimate = new("DistMat3D", values = cor_data[2,],
                              ncol = dim(x$nri)[1], nlyr = 1) 
           )
  attr(final,"is.predictor.nri") <- nri_response
  attr(final,"function") <- "cor.test"
  res@multivariate <- final
  return(res)
}
)
