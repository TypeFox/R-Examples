
gpd.fit <-  function(x,method)
{
    if(!(method %in% c("amle","combined")))
      stop("Unknown method. Please check the ``method'' argument in the help files.")
    if( min(x) < 0 )
	   stop("There are negative observations. \nAll data must be positive real numbers.")
    n <- length(x)
    fit <- switch(method,'amle'=.amle(x,n, ceiling(.2*n)),
                  'combined' = .combined(x)
                  )
    fit <- as.matrix(fit,ncol=1)
    colnames(fit) <- c("Parameter estimate")
    rownames(fit) <- c("shape","scale")
    return(fit)
}
