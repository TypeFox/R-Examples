#'Computes completeness values of the dataset
#'
#'Computes completeness values for each cell. Currently returns Chao2 index of 
#'species richness.
#'
#'After dividing the extent of the dataset in cells (via the
#'\code{\link{getcellid}} function), the function calculates the Chao2 estimator
#'of species richness. Given the nature of the calculations, a minimum number of
#'records must be present on each cell to properly compute the index. If there
#'are too few records in the cells, the function is unable to finish, and it
#'throws an error.
#'
#'@import sqldf
#'@importFrom stats na.omit
#'@importFrom graphics plot
#'@param indf input data frame containing biodiversity data set
#'@param recs minimum number of records per grid cell required to make the 
#'  calculations. Default is 50. If there are too few records, the function 
#'  throws an error.
#'@return data.frame with the columns \itemize{ \item{"Cell_id"}{id of the cell}
#'  \item{"Sobs"}{Number of Observed species} \item{"Sest"}{Estimated number of 
#'  species} \item{"c"}{Completeness ratio the cell} Plots a graph of Number of 
#'  species vs completeness }
#'@examples \dontrun{
#'bdcomplete(inat)
#'}
#'@seealso \code{\link{getcellid}}
#'@export
bdcomplete <- function(indf,recs=50){
  dat1=sqldf("select Scientific_name, Date_collected, Cell_id from indf group by Scientific_name, Date_collected, Cell_id")
  dat2=sqldf("select cell_id,count(*) as cell_ct from dat1 group by cell_id")
  dat3=sqldf(paste("select * from dat2 where cell_ct > ",recs))
  dat1=na.omit(dat1)
  dat2=na.omit(dat2)
  dat3=na.omit(dat3)
  retmat=NULL
  if(dim(dat3)[1]<1){
    stop("Too few data records to compute completeness")
    #return(NULL)
  }
  for (i in 1:dim(dat3)[1]){
    Cell_id=dat3$Cell_id[i]
    cset = dat1[which(dat1$Cell_id==dat3$Cell_id[i]),]
    csum = sqldf("select Scientific_name, count(*) as sct from cset group by Scientific_name")
    Q1=as.numeric(sqldf("select count(*) from csum where sct = 1 "))
    Q2=sqldf("select count(*) from csum where sct = 2 ")
    m=sqldf("select count(*) from ( select * from cset group by Date_collected )")
    Sobs = as.numeric(sqldf("select count(*) from ( select * from cset group by Scientific_name)"))
    Sest = as.numeric(Sobs + (((m-1)/m) * ((Q1 *(Q1 -1))/(2 * (Q2+1)))))
    c = Sobs / Sest
    retmat=rbind(retmat,c(Cell_id,Sobs,Sest,c))
  }
  retmat=as.data.frame(retmat)
  names(retmat)=c("Cell_id","Sobs","Sest","c")
  plot(retmat$Sobs,retmat$c,main="Completeness vs number of species",xlab="Number of species",ylab="Completeness")
  return(retmat)
}