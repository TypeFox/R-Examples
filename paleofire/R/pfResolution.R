#' Calculates age resolution indicators for charcoal records
#' 
#' Calculates age resolution indicators for charcoal records selected using
#' \code{\link{pfSiteSel}} or \code{\link{pfInteractive}} functions.
#' 
#' 
#' @param ID An object of the class "pfSiteSel"
#' @param AgeLim Numeric, defines age limits for age resolution calculations
#' (e.g. AgeLim=c(-50,6000))
#' @return
#' 
#' \item{data.frame}{A data frame with the following informations: ID_SITE,
#' SITE_NAME, Median Resolution of the record, Mean Resolution and Standard
#' deviation}
#' @author O. Blarquez
#' @examples
#' 
#' ID=pfSiteSel(lat>40, lat<90, long>-100, long<=-50)
#' Res=pfResolution(ID,AgeLim=c(-50,8000))
#' head(Res)
#' 
pfResolution=function(ID,AgeLim=NULL){
  
  # Temporal resolution of records
  
  paleofiredata=NULL; rm(paleofiredata)
  
  data(paleofiredata, envir = environment())
  IDs=ID # Save for output
  
  ID=ID$id_site
  # Use only paleofiredata corresponding to ID
  if (is.null(AgeLim)){
  paleofiredata=paleofiredata[paleofiredata[,1] %in% ID,]
  } else {
    paleofiredata=paleofiredata[paleofiredata[,1] %in% ID,]
    paleofiredata=paleofiredata[paleofiredata$EST_AGE>AgeLim[1] & paleofiredata$EST_AGE<AgeLim[2],]
  }
  
  meanres=c()
  medianres=c()
  sdres=c()
  
  for (i in 1:length(ID)){
    meanres[i]=mean(diff(paleofiredata[paleofiredata$ID_SITE %in% ID[i],]$EST_AGE))
    medianres[i]=median(diff(paleofiredata[paleofiredata$ID_SITE %in% ID[i],]$EST_AGE))
    sdres[i]=sd(diff(paleofiredata[paleofiredata$ID_SITE %in% ID[i],]$EST_AGE))
  }
  
  res=data.frame(ID_SITE=as.numeric(IDs$id_site),
                 SITE_NAME=IDs$site_name,
                 MeanRes=as.numeric(meanres),
                 MedianRes=as.numeric(medianres),
                 SdRes=as.numeric(sdres))
  return(res)
  
}
