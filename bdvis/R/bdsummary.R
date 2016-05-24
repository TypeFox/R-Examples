#'Provides summary of biodiversity data
#'
#'Calculates some general indicators of the volume, spatial, temporal and
#'taxonomic aspects of the provided data set.
#'
#'The function returns information on the volume of the data set (number of
#'records), temporal coverage (minimum and maximum dates), taxonomic coverage
#'(brief breakdown of the records by taxonomic levels) and spatial coverage
#'(coordinates of the edges of the bounding box containing all records and
#'division of covered area in degree cells) of the records.
#'
#'@import sqldf
#'@param indf input data frame containing biodiversity data set
#'@export
#'@family Data preparation functions
#'@examples \dontrun{
#'require(rinat)
#'inat<-get_inat_obs_project("reptileindia") 
#'bdsummary(inat)
#'}
bdsummary <- function(indf){
  if(is.null(indf) || is.na(indf) || dim(indf)[1]==0){
    cat("No data in dataset \n")
    return()
  }
  cellcover=sqldf("select count(*) from (select * from indf group by cell_id)")
  names(indf)=gsub("\\.","_",names(indf))
  cat(paste("\n Total no of records =",dim(indf)[1],"\n"))
  cat(paste("\n Temporal coverage... \n"))
  if("Date_collected" %in% colnames(indf)){
    cat(paste(" Date range of the records from ",
              range(as.Date(indf$Date_collected),na.rm=T)[1]," to ",
              range(as.Date(indf$Date_collected),na.rm=T)[2],"\n"))
  }
  cat(paste("\n Taxonomic coverage... \n"))
  cat(paste(" No of Families : ",length(unique(indf$Family)),"\n"))
  cat(paste(" No of Genus : ",length(unique(indf$Genus)),"\n"))
  cat(paste(" No of Species : ",length(unique(indf$Scientific_name)),"\n"))
  cat(paste("\n Spatial coverage ... \n"))
  cat(paste(" Bounding box of records ",min(indf$Latitude,na.rm=T),",",min(indf$Longitude,na.rm=T),
            " - ",max(indf$Latitude,na.rm=T),",",max(indf$Longitude,na.rm=T),"\n"))
  if (!("cell_id" %in% colnames(indf))){
    indf=getcellid(indf)
  }
  cat(paste(" Degree celles covered : ",cellcover,"\n"))
  latrng=ceiling(max(indf$Latitude,na.rm=T))-ceiling(min(indf$Latitude,na.rm=T))
  longrng=ceiling(max(indf$Longitude,na.rm=T))-ceiling(min(indf$Longitude,na.rm=T))
  cat(paste(" % degree cells covered : ",(cellcover/(latrng * longrng))*100,"\n"))
}
