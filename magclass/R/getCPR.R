getCPR <- function(x) {
  if(!is.magpie(x)) {
    if(x==0.5) {
      region.code <- magclassdata$half_deg$region
      cpr <- rep(0,length(levels(region.code)))
      names(cpr) <- levels(region.code)
      for(region in names(cpr)) {
       cpr[region] <- length(grep(region,region.code))  
      }
    } else {
      stop(paste("No cells-per-region information available for resolution",x))
    }
  } else {
    region_names <- getRegions(x)
    cpr <- rep(0,length(region_names))
    names(cpr) <- region_names
    for(region in region_names) {
      cpr[region] <- length(grep(region,dimnames(x)[[1]]))  
    }
  }
  return(cpr)
}