#' Distribution graphs
#' 
#' Build plots dispalying distribution of biodiversity records among 
#' user-defined features.
#' 
#' The main use of this function is to create record histograms according to
#' different features of the data set. For example, one might want to see the 
#' evolution of records by year, or by species. This function enables easy 
#' access to such plots.
#' 
#' @import sqldf
#' @importFrom graphics hist plot
#' @param indf input data frame containing biodiversity data set
#' @param ptype Feature to represent. Accepted values are "species", "cell" and 
#'   "efforts" (year)
#' @param ... any additional parameters for the \code{\link{plot}} function.
#' @examples \dontrun{
#'  distrigraph(inat,ptype="cell",col="tomato")
#'  distrigraph(inat,ptype="species",ylab="Species")
#'  distrigraph(inat,ptype="efforts",col="red")
#'  distrigraph(inat,ptype="efforts",col="red",type="s")
#' }
#' @family Data preparation functions
#' @export
distrigraph <- function(indf,ptype=NA,...){
  custgraph='col="red"'
  if(!is.na(ptype)){
    switch(ptype,
           cell={
             mat=sqldf("select Cell_id, count(*) as Records from indf group by Cell_id")
             hist(mat$Records,main="Distribution of Records per cell",xlab="Records",...)
           },
           species={
             mat=sqldf("select Scientific_name, count(*) as Records from indf group by Scientific_name")
             hist(mat$Records,main="Distribution of Records per Species",xlab="Records",...)
           },
           efforts={
             Year = as.numeric(strftime(as.Date(indf$Date_collected,na.rm=T), format = "%Y"))
             indf=cbind(indf,Year)
             mat=sqldf("select Year, count(*) as Records from indf group by Year order by Year")
             plot(mat,main="Distribution of collection effforts over time",...)
           },
           stop("Not a valid option. See ?distrigraph for currently accepted values")
    )
  } else {
    stop("Must indicate a value for ptype.")
  }
}