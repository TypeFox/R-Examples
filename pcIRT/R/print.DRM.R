#' @rdname drm
#'
#' @export
#'
#' @method print DRM
#'
#' @param x object of class \code{DRM}

print.DRM <-
  function(x, ...){
    
    parall <- rbind(cbind("item estimates"=x$itempar, "SE"=x$itempar_se, "SE low"= x$itempar_se_low, "SE up"=x$itempar_se_up),"distrpar"=c(x$distrpar, x$distrpar_se, x$distrpar_se_low, x$distrpar_se_up) )
    
    cat("Parameter estimates: \n")
    print(parall)  
    
  }
