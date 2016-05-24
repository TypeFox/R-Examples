# Wrapper function for outlier plots
outlierPlot <- function(y, L, mode="qq", ...)
{
   if ( mode=="residual" & L$method != "Method II" )
      stop("resiudal plot is only available for Method II")
   
   switch(mode,
      qq = qqFitPlot(y, L, ...),
      residual = plotMethodII(y, L, ...),
      stop("Plot mode not recognized")
      )

}


