#' @title Converting a booklet allocation table into a incidence matrix
#' @export make.incidenz
#' @description This function converts a booklet allocation table (like in \code{\link{cogBOOKLET}}) into a incidence matrix used in the function \code{\link{pers}}.
#' @details It is assumed that there is an equal replicate factor for each item used, when constructing the bookletdesign - so every items occures with the same frequency over all booklets of the entire set of booklets.
#'
#' @param tab a booklet allocation table as a \code{data.frame}. The first column is assumed to contain the item names as a character vector (not a factor!) the other columns must be integer vectors containing the information in which booklet(s) the respective item is allocated.
#' @param bookid a integer vector with the same length as the number of persons in the response data giving the information which booklet was assigned to each person.  
#' @param item_order optional a character vector with the item names in the order of the items in the response data (from first to last column in the response data). By default it is assumend that the item order in the booklet allocation table is already the same as in the response data.
#' @param info logical default: \code{info=FALSE} to return just the incidence matrix. If set to \code{info=TRUE} more detailed information about the booklet design ist returned. 
#' @return an incidence matrix as an object of class "matrix" with 0,1 coding or a "list" with detailed information.

#' @examples 
#' #########################
#' data(cog);data(cogBOOKLET) # loading reponse and allocation data
#' table(cog$BOOKID)# show n persons per booklet
#' names(table(c(as.matrix(cogBOOKLET[,2:5])))) # show booklets in allocation data
#' d<-(cog[cog$BOOKID!=14,]) # skip persons which got booklet No.14.
#' inc<-make.incidenz(tab=cogBOOKLET, bookid=d$BOOKID) # make just the incidence matrix
#' inc  
#' make.incidenz(tab=cogBOOKLET, bookid=d$BOOKID, info=TRUE) # get some info too
#' # in this case not necessary but just to show
#' # using the (item) names in cog to secure the item order in incidence matrix:
#' make.incidenz(tab=cogBOOKLET, bookid=d$BOOKID, item_order=names(cog)[4:34])  
#' #######################

make.incidenz <- function(tab, bookid, item_order=NULL, info=FALSE) {
  n_items<-dim(tab)[1]
  n_repli<-dim(tab)[2]-1
  n_booklets<-length(unique(bookid))
  n_persons<-length(bookid)
  alloc<-tab[ ,2:(n_repli+1)]
  
  boolets_alloc<-sort(unique(c(as.matrix(alloc))))
  boolets_bookid<-sort(unique(bookid))
  
  stopifnot(n_booklets==max(bookid))
  if( (length(boolets_alloc) != length(boolets_bookid)) & !all(boolets_bookid%in%boolets_alloc)   ){
    cat("number and naming of boklets in 'tab' and 'bookid' dont match !","\n")
    stop
  }
  
  if (length(item_order)!=0){
    stopifnot(length(item_order)==n_items)
  }
  if (length(item_order)!=0){ # optionales sortieren von tab nach item_order
    rownames(tab)<-tab[,1]
    tab<-tab[item_order,] 
  }
  # matrix items nach booklets -------
  item_booklet <- matrix(NA,nrow=n_items,ncol=n_booklets)
  rownames(item_booklet)<-tab[,1]
  colnames(item_booklet)<-paste("Booklet",formatC(1:n_booklets, width = nchar(paste(n_booklets)), format = "d", flag = "0"),sep=" ")
  for (i in 1:n_booklets){
    item_booklet[,i]<- rowSums(alloc==i)  
  } 
  # matrix pesonen nach items (nach booklets) -------  
  incidenz <- matrix(NA,nrow=n_persons,ncol=n_items)
  colnames(incidenz)<-tab[,1]
  for (i in 1:n_persons){
    incidenz[i,] <-      item_booklet[  ,bookid[i] ]
  } 
  if(info==TRUE){
    result<-list(incidenz=incidenz,item_booklet=item_booklet,replicate_factor=n_repli,n_booklets=n_booklets, table_booklets=table(as.matrix(alloc))  )
  }
  if(info==FALSE){
    result<-incidenz } 
  return(result)
}
