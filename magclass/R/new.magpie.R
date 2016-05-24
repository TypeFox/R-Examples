new.magpie <- function(cells_and_regions="GLO",years=NULL,names=NULL,fill=NA,sort=FALSE,sets=NULL) {
  ncells <- length(cells_and_regions)
  nyears <- ifelse(is.null(years),1,length(years))
  ndata  <- ifelse(is.null(names),1,length(names))
  if(all(!grepl("\\.",cells_and_regions)) & !is.null(cells_and_regions)) {
    if(all(is.numeric(cells_and_regions))) {
      cells_and_regions <- paste("GLO",cells_and_regions,sep=".")
    } else {
      #cells_and_regions <- paste(cells_and_regions,1:ncells,sep=".")
    }
  }
  object<-new("magpie",array(fill,dim=c(ncells,nyears,ndata)))
  getCells(object) <- cells_and_regions
  getYears(object) <- years
  getNames(object) <- names
  if(sort) object <- magpiesort(object)
  object <- clean_magpie(object,"sets")
  if(!is.null(sets)) getSets(object) <- sets
  return(object)
}