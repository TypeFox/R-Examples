data.frc <- function(data.in,method=c("crost","crost.ma","tsb","sexsm","imapa","auto"),...){
# Wrapper to forecasts data.frames with a single call
# 
# Inputs:
#   data.in     Data frame with time series. This can also be a matrix or array with each 
#               column being a different time series. 
#   method      Which method to use for forecasting: 
#                 "crost", "crost.ma", "tsb", "sexsm", "imapa", "auto" 
#               "auto" use PKa classification to select between Croston, SBA and SES.
#   ...         Additional inputs to pass to forecasting functions. See individual function 
#               documentation for options. 
#
# Outputs:
#   frc.out     Data frame containing forecasts for all time series.
#   out         List with detailed output per series. To access individual outputs of the list
#               use: sapply(out, get, x="element"), where "element" could be for example "frc.in". 
#
# Nikolaos Kourentzes, 2015 <nikolaos@kourentzes.com>

  # Get defaults
  method <- method[1]
  
  # Get number of columns and class of data.in
  # k <- ncol(data.in)
  # data.cl <- class(data.in)
  data.names <- colnames(data.in)
  
  # Check is selected model is implemented
  allow.method <- c("crost","crost.ma","tsb","sexsm","imapa","auto")
  if (!(method %in% allow.method)){
    stop(paste(method,"not a permitted forecasting method."))
  }
  
  # Fit models and produce forecasts
  switch(method,
         "crost" = {out <- apply(data.in,2,crost,...)},
         "crost.ma" = {out <- apply(data.in,2,crost.ma,...)},
         "tsb" = {out <- apply(data.in,2,tsb,...)},
         "sexsm" = {out <- apply(data.in,2,sexsm,...)},
         "imapa" = {out <- apply(data.in,2,imapa,...)},
         "auto" = # Select between Croston, SBA and SES
           {
             out <- vector("list",sum(cls$summary))
             # Distribute elipse (...) inputs to each method
             dots <- list(...)
             valid.crost <- c("h","w","init","nop","cost","init.opt","outplot","na.rm")
             dots.crost <- dots[names(dots) %in% valid.crost]
             valid.sexsm <- c("h","w","init","cost","init.opt","outplot","na.rm")
             dots.sexsm <- dots[names(dots) %in% valid.sexsm]
             # Perform classification
             cls <- idclass(data.in,type="PKa",outplot="none")
             # Forecast
             if (length(cls$idx.croston)>=1){
               # out.croston <- apply(data.in[,cls$idx.croston],2,do.call(crost,c(list(type="croston"),dots.crost)))
               out.croston <- do.call(function(...){apply(data.in[,cls$idx.croston],2,crost,...)}, c(list(type="crost"),dots.crost))
               out[cls$idx.croston] <- out.croston
             } else {
               out.croston <- NULout.croston <- do.call(function(...){apply(data.in[,cls$idx.croston],2,crost,...)}, c(list(type="sba"),dots.crost))
             }
             if (length(cls$idx.sba)>=1){
               # out.sba <- apply(data.in[,cls$idx.sba],2,crost,type="sba",dots.crost) 
               out.sba <- do.call(function(...){apply(data.in[,cls$idx.sba],2,crost,...)}, c(list(type="sba"),dots.crost))
               out[cls$idx.sba] <- out.sba
             } else {
               out.sba <- NULL
             }
             if (length(cls$idx.ses)>=1){
               # out.ses <- apply(data.in[,cls$idx.ses],2,sexsm,dots.sexsm)
               out.ses <- do.call(function(...){apply(data.in[,cls$idx.ses],2,sexsm,...)}, dots.sexsm)
               out[cls$idx.ses] <- out.ses
             } else {
               out.ses <- NULL
             }
             names(out) <- data.names
           }
         )
  
  # Now get from res all forecasts and return them to a data.frame
  frc.out <- sapply(out, get, x="frc.out")
  colnames(frc.out) <- data.names
  frc.out <- data.frame(frc.out)
  
  return(list(frc.out = frc.out, out = out))

}