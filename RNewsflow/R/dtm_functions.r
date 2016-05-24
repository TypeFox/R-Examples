dtmToSparseMatrix <- function(dtm){
  sm = Matrix::spMatrix(nrow(dtm), ncol(dtm), dtm$i, dtm$j, dtm$v)
  rownames(sm) = rownames(dtm)
  colnames(sm) = colnames(dtm)
  sm
}

confirm.dtm.meta <- function(meta, id.var, date.var){
  if(!id.var %in% colnames(meta)) stop(sprintf('Meta data.frame should contain a column that matches the id.var parameter (currently set to "%s")', id.var))
  if(!date.var %in% colnames(meta)) stop(sprintf('Meta data.frame should contain a column that matches the id.var parameter (currently set to "%s")', date.var))
}

match.dtm.meta <- function(dtm, meta, id.var){
  if(mean(rownames(dtm) %in% meta[,id.var]) < 1) stop('Not all documents in DTM match with a document in the meta data.frame')
  meta[match(rownames(dtm), meta[,id.var]),]
}


#' Calculate statistics for term occurence across days
#'
#' @param dtm A document-term matrix in the tm \link[tm]{DocumentTermMatrix} class or a TsparseMatrix from the Matrix class (\link[Matrix]{spMatrix}) 
#' @param meta A data.frame where rows are documents and columns are document meta information. 
#' Should contain 2 columns: the document name/id and date. 
#' The name/id column should match the rownames (i.e. document names) of the DTM, and its label is specified in the `id.var` argument. 
#' The date column should be intepretable with \link[base]{as.POSIXct}, and its label is specified in the `date.var` argument.            
#' @param id.var The label for the document name/id column in the `meta` data.frame. Default is "document_id"
#' @param date.var The label for the document date column in the `meta` data.frame . default is "date"
#'
#' @return A data.frame with statistics for each term.
#' \itemize{
#'  \item{freq:}{ The number of times a term occurred}
#'  \item{doc.freq:}{ The number of documents in which a term occured}
#'  \item{days.n:}{ The number of days on which a term occured}
#'  \item{days.pct:}{ The percentage of days on which a term occured}
#'  \item{days.entropy:}{ The entropy of the distribution of term frequency across days}
#'  \item{days.entropy.norm:}{ The normalized days.entropy, where 1 is a discrete uniform distribution}
#' }
#' @export
#'
#' @examples
#' data(dtm)
#' data(meta)
#' 
#' tdd = term.day.dist(dtm, meta)
#' head(tdd)
#' tail(tdd)
term.day.dist <- function(dtm, meta, id.var='document_id', date.var='date'){
  confirm.dtm.meta(meta, id.var, date.var)
  meta = match.dtm.meta(dtm, meta, id.var)
  
  if('DocumentTermMatrix' %in% class(dtm)) dtm = dtmToSparseMatrix(dtm)
  
  cs = Matrix::colSums(dtm)
  if(sum(cs == 0) > 0) {
    message("dtm contains empty columns/terms. These will be ignored (and won't appear in the output)")
    dtm = dtm[,slam::col_sums(dtm) > 0]
  } 

  document.date = as.Date(meta[,date.var])
  dateseq = seq.Date(min(document.date), max(document.date), by='days')
  i = document.date[dtm@i+1]
  i = match(i, dateseq)
  m = Matrix::spMatrix(length(dateseq), ncol(dtm), i, dtm@j+1, dtm@x)
  
  m = methods::as(m, 'dgCMatrix')
  days.entropy = columnEntropy(m)
  days.n = Matrix::colSums(m>0)

  d = data.frame(term=colnames(dtm),
                 freq = Matrix::colSums(dtm),
                 doc.freq = Matrix::colSums(dtm > 0),
                 days.n = days.n, 
                 days.pct = days.n / length(dateseq),
                 days.entropy = days.entropy, 
                 days.entropy.norm=days.entropy / length(dateseq))
  d$term = as.character(d$term)
  rownames(d) = NULL
  d
}
