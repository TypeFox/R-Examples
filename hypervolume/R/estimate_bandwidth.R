estimate_bandwidth <- function(data,method="silverman")
{
  # remove missing cases
  data <- na.omit(data)
  npoints = nrow(data)
  ndim <- ncol(data)
  
  if (method=="silverman")
  {
    stdev = apply(data,2,sd)
    
    normal_silverman_bw = 1.06 * stdev * npoints ^ (-1/5)
    
    # scale to box width
    rect_bw = 2*normal_silverman_bw
    
    # remove names
    names(rect_bw) <- NULL
    
    return(rect_bw)    
  }
  else
  {
    if (ndim <= 6)
    {
     if (method=="plug-in")
     {
       bw_plugin <- ks::Hpi.diag(data, nstage=2, pilot="samse",pre="scale")
       
       # convert estimated variances to standard deviations and then box-widths
       bw_plugin_scaled <- 2*sqrt(diag(bw_plugin))
       
       return(bw_plugin_scaled)
     }
     else if (method=="cross-validation")
     {
       bw_crossvalidation <- ks::Hscv.diag(data, pilot="samse", pre="scale")
       
       # convert estimated variances to standard deviations and then box-widths
       bw_crossvalidation_scaled <- 2*sqrt(diag(bw_crossvalidation))
       
       return(bw_crossvalidation_scaled)
     }
     else
     {
       stop("Argument 'method' not recognized")
     }
    }
    else
    {
      stop("Non-Silverman bandwidth estimators not available for n>6 dimensions.")
    }
  }
}