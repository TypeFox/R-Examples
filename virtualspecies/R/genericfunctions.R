#' @export
#' @method print virtualspecies
print.virtualspecies <- function(x, ...)
{
  cat(paste("Virtual species generated from", 
            length(x$details$variables), 
            "variables:\n",
            paste(x$details$variables, collapse = ", ")))
  cat("\n\n- Approach used:")
  if(x$approach == "response")
  {
    cat(" Responses to each variable")

    cat("\n- Response functions:")
    sapply(x$details$variables, FUN = function(y)
    {
      cat("\n   .", y, 
          "  [min=", x$details$parameters[[y]]$min, 
          "; max=", x$details$parameters[[y]]$max,
          "] : ", 
          x$details$parameters[[y]]$fun,
          "   (", 
          paste(names(x$details$parameters[[y]]$args), 
                x$details$parameters[[y]]$args, sep = '=', collapse = "; "),
          ")", sep = "")
    })
    
    if (x$details$rescale.each.response)
    {
      cat("\n- Each response function was rescaled between 0 and 1")
    } else
    {
      cat("\n- Response functions were not rescaled between 0 and 1")
    }
    
    cat("\n- Environmental suitability formula = ", x$details$formula, sep = "")
    
    if (x$details$rescale)
    {
      cat("\n- Environmental suitability was rescaled between 0 and 1")
    } else
    {
      cat("\n- Environmental suitability was not rescaled between 0 and 1")
    }
    
  } else if(x$approach == "pca")
  {
    cat(" Response to axes of a PCA")
    cat("\n- Axes: ", 
        paste(x$details$axes, collapse = " & "),
        "; ", round(sum(x$details$pca$eig[x$details$axes])/sum(x$details$pca$eig) * 100, 2),
        "% explained by these axes")
    cat("\n- Responses to axes:")
    sapply(c(1, 2), function(y)
    {
      cat("\n   .Axis ", x$details$axes[y],
          "  [min=", round(min(x$details$pca$li[, x$details$axes[y]]), 2),
          "; max=", round(max(x$details$pca$li[, x$details$axes[y]]), 2),
          "] : dnorm (mean=", x$details$means[y], "; sd=", x$details$sds[y],
          ")", sep = "")
    })
    if (x$details$rescale)
    {
      cat("\n- Environmental suitability was rescaled between 0 and 1")
    } else
    {
      cat("\n- Environmental suitability was not rescaled between 0 and 1")
    }
  }
  if(!is.null(x$PA.conversion))
  {
    cat("\n- Converted into presence-absence:")
    cat("\n   .Method =", x$PA.conversion["conversion.method"])
    if(x$PA.conversion["conversion.method"] == "probability")
    {
      cat("\n   .alpha (slope)           =", x$PA.conversion["alpha"])
      cat("\n   .beta  (inflexion point) =", x$PA.conversion["beta"])
      cat("\n   .species prevalence      =", x$PA.conversion["species.prevalence"])
    } else if(x$PA.conversion["conversion.method"] == "threshold")
    {
      cat("\n   .threshold           =", x$PA.conversion["cutoff"])
      cat("\n   .species prevalence  =", x$PA.conversion["species.prevalence"])
    }
  }
  if(!is.null(x$occupied.area))
  {
    if(!is.null(x$geographical.limit))
    {
      cat("\n- Distribution bias introduced:")
      cat("\n   .method used :", x$geographical.limit$method)
      if(x$geographical.limit$method %in% c("country", "region", "continent"))
      {
        cat("\n   .area(s)     :", x$geographical.limit$area)
      } else if(x$geographical.limit$method == "extent")
      {
        cat("\n   .extent      : [Xmin; Xmax] = [", 
            x$geographical.limit$extent@xmin, "; ",
            x$geographical.limit$extent@xmax, "] - [Ymin; Ymax] = [",
            x$geographical.limit$extent@ymin, "; ",
            x$geographical.limit$extent@ymax, "]", sep = "")
      } else if(x$geographical.limit$method == "polygon")
      {
        cat("\n   .polygon    : Object of class ", class(x$geographical.limit$area), sep = "")
      }
    }
  }
}

#' @export
#' @method str virtualspecies
str.virtualspecies <- function(object, ...)
{
  args <- list(...)
  if(is.null(args$max.level))
  {
    args$max.level <- 2
  }
  NextMethod("str", object = object, max.level = args$max.level)
}

#' @export
#' @method plot virtualspecies
plot.virtualspecies <- function(x, ...)
{
  y <- raster::stack(x$suitab.raster)
  names(y) <- "Suitability.raster"
  if(!is.null(x$pa.raster))
  {
    y <- stack(y,
               x$pa.raster)
    names(y)[[nlayers(y)]] <- "Presence.absence.raster"
  }
  if(!is.null(x$occupied.area))
  {
    y <- stack(y,
               x$occupied.area)
    names(y)[[nlayers(y)]] <- "Occupied.area.raster"
  }
  x <- y
  plot(x, ...)
}