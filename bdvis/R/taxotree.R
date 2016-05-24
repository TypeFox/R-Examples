#'Treemap based on taxonomic hierarchy of records
#'
#'Draws a treemap (\url{https://en.wikipedia.org/wiki/Treemapping}) based on the
#'taxonomic information of the records.
#'
#'This function builds a treemap of the taxonomic information present in the 
#'data set. It represents this information at two levels (with the arguments 
#'sum1 and sum2). The first level (sum1) will be represented with cell sizes and
#'is a reflection of the number of records in that group. If, for example, 
#'"Family" is selected as value for sum1, the size of the cells in the treemap 
#'will be directly proportional to the number of records for that taxonomic 
#'family. The second level (sum2) will be represented by color and is a 
#'reflection of the number of sub-groups in a particular cell. If, for example, 
#'"Genus" is selected as value for sum2, the color of the cell will depend on 
#'the number of different genera for that particular cell.
#'
#'@import sqldf
#'@import treemap
#'@importFrom utils tail
#'@param indf input data frame containing biodiversity data set
#'@param n maximum number of rectangles to be plotted in the treemap. Default is
#'  30
#'@param title title for the tree. Default is "Records per <sum1>"
#'@param legend legend title. Default is "Number of <sum2>"
#'@param sum1 Taxonomic level whose density will be represented with different 
#'  cell sizes
#'@param sum2 Taxonomic level whose density will be represented with a color 
#'  gradient
#'@references Otegui, J., Arino, A. H., Encinas, M. A., & Pando, F. (2013). 
#'  Assessing the Primary Data Hosted by the Spanish Node of the Global 
#'  Biodiversity Information Facility (GBIF). PLoS ONE, 8(1), e55144. 
#'  doi:10.1371/journal.pone.0055144
#'@family Taxonomic visualizations
#'@export
#'@examples \dontrun{
#'taxotree(inat)
#'}
taxotree <- function(indf,n=30,title=NA,legend=NA,sum1="Family",sum2="Genus"){
  if(!is.element(sum1,names(indf))){
    cat(paste("Field",sum1," not found in dataset \n"))
    return()
  }
  if(!is.element(sum2,names(indf))){
    cat(paste("Field",sum2," not found in dataset \n"))
    return()
  }
  if (!is.na(title)) {
    title2 <- title
  } else {
    title2 <- paste("Records per",sum1)
  }
  if (!is.na(legend)) {
    legend2 <- legend
  } else {
    legend2 <- paste("Number of",sum2)
  }
  sql1=paste("select",sum1,",",sum2,",count(*) as Recs from indf group by ",sum1,",",sum2 )
  tab1=(sqldf(sql1))
  sql2=paste("select",sum1,",count(*) as Gnum, sum(Recs) as Rec1 from tab1 group by ",sum1,"order by Rec1")
  tab2=sqldf(sql2)
  tab2=tail(tab2,n)
  treemap(tab2, index=c(sum1), vSize="Rec1", vColor="Gnum", type="value",
          title=title2, title.legend=legend2)
}