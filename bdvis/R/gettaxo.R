#'Get higher taxonomy data
#'
#'Retrieve higher taxonomy information (like Family and Order) for each record 
#'from the "Encyclopedia of Life" web API.
#'
#'This function makes use of certain functions in the \code{\link{taxize}}
#'package. It scans and retrieves the taxonomic hierarchy for each scientific
#'name (or just genus name) in the data set. When new data are retrieved, they
#'are stored in a local sqlite database, taxo.db, for faster further access.
#'
#'@import plyr
#'@import taxize
#'@import sqldf
#'@importFrom utils setTxtProgressBar txtProgressBar
#'@param indf input data frame containing biodiversity data set
#'@param genus If TRUE, use only genus level data to get taxanomy
#'@param verbose If TRUE, displays each name string for which the higher 
#'  taxonomy is sought
#'@return indf with added / updated columns \itemize{ \item{"Kingdom"}{Kingdom 
#'  of the Scientific name} \item{"Phylum"}{Phylum of the Scientific name} 
#'  \item{"Order_"}{Order of the Scientific name} \item{"Family"}{Family of the 
#'  Scientific name} \item{"Genus"}{Genus of the Scientific name} } and also 
#'  saves a local copy of taxanomy downloaded for future use in taxo.bd sqlite 
#'  file
#'@examples \dontrun{
#'inat=gettaxo(inat)
#'}
#'@family Data preparation functions
#'@export
gettaxo <- function(indf,genus=FALSE,verbose=FALSE){
  names(indf)[names(indf)=="Order_"]<-"Order"
  indfu=sqldf("select Scientific_name from indf group by Scientific_name order by Scientific_name")
  if(file.exists("taxo.db")){
    mytaxo=sqldf("select * from taxo", dbname = "taxo.db")
    indfu=sqldf("select indfu.*, mytaxo.* from indfu left outer join mytaxo on  indfu.Scientific_name = mytaxo.Scientific_name")
    indfu=indfu[,-c(2)]
    cat("Local database taxo.db is present\n")
  } else {
    indfu$Kingdom[1]=NA
    indfu$Phylum[1]=NA
    indfu$Class[1]=NA
    indfu$Order[1]=NA
    indfu$Family[1]=NA
    indfu$Genus[1]=NA
    cat("Local database taxo.db is absent\n")
  }
  indfu1=indfu[which(!is.na(indfu$Kingdom)),]
  indfu=indfu[which(is.na(indfu$Kingdom)),]
  if(dim(indfu)[1]>0){
    cat("Processing names not present in local database ...\n")
    sciname=""
    dat1<-NULL
    for(i in 1:dim(indfu)[1]){
      if(dim(indfu)[1]<=1){pbar=2} else {pbar=dim(indfu)[1]}
      pb   <- txtProgressBar(1, pbar, style=3)
      if(genus){
        Scientific_name=strsplit(indfu$Scientific_name[i]," ")[[1]][1]
        if(is.na(Scientific_name)){Scientific_name=""}
      } else {
        Scientific_name=indfu$Scientific_name[i]
      }
      if (!Scientific_name==sciname){
        dat <- classification(get_uid(Scientific_name,verbose = verbose))
        if(!is.na(dat)){
          dat1<-NULL
          dat1 <- ldply(dat, function(x) x[x$rank %in% c("kingdom","phylum","class","order","family","genus","species"), "name"])
          dat2 <- ldply(dat, function(x) x[x$rank %in% c("kingdom","phylum","class","order","family","genus","species"), "rank"])
        }
      }
      if(!is.null(dat1)) {
        names(dat1) <- dat2
        if(!is.null(dat1$kingdom)) indfu$Kingdom[i] <- dat1$kingdom
        if(!is.null(dat1$phylum)) indfu$Phylum[i] <- dat1$phylum
        if(!is.null(dat1$class)) indfu$Class[i] <- dat1$class
        if(!is.null(dat1$order)) indfu$Order[i] <- dat1$order
        if(!is.null(dat1$family)) indfu$Family[i] <- dat1$family
        if(!is.null(dat1$genus)) indfu$Genus[i] <- dat1$genus
      }
      sciname <- Scientific_name
      if(!verbose) {setTxtProgressBar(pb, i)}
    }
  }
  if(!verbose) {setTxtProgressBar(pb, pbar)}
  indfu=rbind(indfu,indfu1)
  indf=sqldf("select * from indf, indfu where indf.Scientific_name = indfu.Scientific_name")
  if(file.exists("taxo.db")){
    cat("\nUpdating local database with new names\n")
    sqldf("insert into taxo select * from indfu where Scientific_name not in (select Scientific_name from taxo)",dbname="taxo.db")
  } else {
    sqldf("attach 'taxo.db' as new")
    sqldf("create table taxo as select * from indfu", dbname="taxo.db")
    cat("\nCreating local database taxo.db ...\n")
  } 
  indf <- indf[, !duplicated(colnames(indf))]
  names(indf)[names(indf)=="Order"]<-"Order_"
  return(indf)
}