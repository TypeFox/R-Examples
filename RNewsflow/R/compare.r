## FUNCTIONS FOR COMPARING DOCUMENTS


cosineSimilarity <- function(m1, m2=NULL){
  norm = sqrt(Matrix::colSums(m1^2))
  m1@x = m1@x / norm[m1@j+1]  
  if(!is.null(m2)){
    norm = sqrt(Matrix::colSums(m2^2))
    m2@x = m2@x / norm[m2@j+1]
    cp = Matrix::crossprod(m1,m2) 
  } else cp = Matrix::crossprod(m1)
  cp
}

termOverlap <- function(m1, m2=m1){
  m2@x[Matrix::which(m2@x > 0)] = 1
  Matrix::crossprod(m1,m2)
}

termOverlap_pct <- function(m1, m2=m1, reverse=FALSE){
  totalterms = if(!reverse) Matrix::colSums(methods::as(m1, 'dgCMatrix')) else Matrix::colSums(methods::as(m2, 'dgCMatrix'))
  m2@x[Matrix::which(m2@x > 0)] = 1
  Matrix::crossprod(m1,m2) / totalterms
}

termIndex <- function(m1, m2=m1){
  m2@x[Matrix::which(m2@x > 0)] = 1
  totalterms = Matrix::colSums(methods::as(m1, 'dgCMatrix'))
  Matrix::crossprod(m1,m2) / totalterms
}

Nth.max <- function(x, N){
  N = min(N, length(x)) 
  -sort(-x, partial=N)[N]
}

filterResults <- function(results, min.similarity, n.topsim){
  if(!is.null(min.similarity)) results@x[Matrix::which(results@x < min.similarity)] = 0
  if(!is.null(n.topsim)) {
    simthres = apply(results, 1, Nth.max, N=n.topsim)
    results@x[Matrix::which(results < simthres)] = 0
  }
  results
}
  
calculate.similarity <- function(m.x, m.y, measure){
  if(measure == 'cosine') results = cosineSimilarity(m.x, m.y)
  if(measure == 'percentage.from') results = termOverlap_pct(m.x, m.y)
  if(measure == 'percentage.to') results = termOverlap_pct(m.x, m.y, reverse = TRUE)
  results
}

reindexTerms <- function(dtm, terms){
  dtm = dtmToSparseMatrix(dtm)
  documents = rownames(dtm)
  dtm = Matrix::spMatrix(nrow(dtm), length(terms), dtm@i+1, match(colnames(dtm)[dtm@j+1], terms), dtm@x)
  dimnames(dtm) = list(documents, terms)
  tm::as.DocumentTermMatrix(dtm, weighting = tm::weightTf)
}

#' Compare the documents in two corpora/dtms
#' 
#' Compare the documents in corpus dtm.x with reference corpus dtm.y. 
#' 
#' The calculation of document similarity is performed using a vector space model approach. 
#' Inner-product based similarity measures are used, such as cosine similarity.
#' It is recommended to weight the DTM beforehand, for instance using Term frequency-inverse document frequency (tf.idf)
#' 
#' @param dtm A document-term matrix in the tm \link[tm]{DocumentTermMatrix} class. It is recommended to weight the DTM beforehand, for instance using \link[tm]{weightTfIdf}.
#' @param dtm.y Optional. If given, documents from dtm will only be compared to the documents in dtm.y
#' @param measure the measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine", for cosine similarity. Also supports assymetrical measures "percentage.from" and "percentage.to" for the percentage of overlapping terms (term scores taken into account). Here "percentage.from" gives the percentage of the document that is compared to the other, whereas "percentage.to" gives the percentage of the document to which is compared.
#' @param min.similarity a threshold for similarity. lower values are deleted. Set to 0 by default.
#' @param n.topsim An alternative or additional sort of threshold for similarity. Only keep the [n.topsim] highest similarity scores for x. Can return more than [n.topsim] similarity scores in the case of duplicate similarities.
#' @param return.zeros If true, all comparison results are returned, including those with zero similarity (rarely usefull and problematic with large data)
#' 
#' @return A data frame with pairs of documents and their similarities. 
#' @export
#' 
#' @import tm
#' 
#' @examples
#' data(dtm)
#' 
#' comp = documents.compare(dtm, min.similarity=0.4)
#' head(comp)
documents.compare <- function(dtm, dtm.y=NULL, measure='cosine', min.similarity=0, n.topsim=NULL, return.zeros=FALSE) {  
  if(!is.null(dtm.y)){
    if(mean(colnames(dtm) == colnames(dtm.y)) < 1){
      ## if colnames do not match, reindex them.
      terms = unique(c(colnames(dtm), colnames(dtm.y)))
      dtm = reindexTerms(dtm, terms)
      dtm.y = reindexTerms(dtm.y, terms)
    }
    m.x = Matrix::t(dtmToSparseMatrix(dtm))
    m.y = Matrix::t(dtmToSparseMatrix(dtm.y))
  } else {
    m.x = m.y = Matrix::t(dtmToSparseMatrix(dtm))
  }
  
  results = calculate.similarity(m.x, m.y, measure)
  results = filterResults(results, min.similarity, n.topsim)
  
  results = methods::as(results, 'dgTMatrix')
  if(return.zeros) {
    results = Matrix::Matrix(Matrix::which(!is.na(results), arr.ind=TRUE))
    results = data.frame(x=colnames(m.x)[results[,1]], y=colnames(m.y)[results[,2]], similarity=as.vector(results))
  } else{
    if(sum(results) == 0) return(NULL)
    results = data.frame(x=colnames(m.x)[results@i+1], y=colnames(m.y)[results@j+1], similarity=results@x)
    results = results[results$similarity > 0 & !is.na(results$similarity),]
  }
  results[!as.character(results$x) == as.character(results$y),]
}

unlistWindow <- function(list_object, i, window){
  indices = i + window
  indices = indices[indices > 0 & indices <= length(list_object)]
  unlist(list_object[indices], use.names=FALSE)
}

getDateIds <- function(date, row_filter=NULL){
  if(is.null(row_filter)) row_filter = rep(TRUE, length(date))
  
  datetime = as.Date(date)
  datetimeseq = seq.Date(min(datetime), max(datetime), by='days')
  
  nonempty = which(datetimeseq %in% unique(datetime))
  nonempty_datetime_ids = plyr::llply(datetimeseq[nonempty], function(dtime) which(datetime == dtime & row_filter))
  datetime_ids = vector("list", length(datetimeseq))
  datetime_ids[nonempty] = nonempty_datetime_ids
  datetime_ids
}

#' Compare the documents in a dtm with a sliding window over time
#' 
#' Given a document-term matrix (DTM) and corresponding document meta data, calculates the document similarities over time using with a sliding window.
#'  
#' The meta data.frame should have a column containing document id's that match the rownames of the DTM (i.e. document names) and should have a column indicating the publication time. 
#' By default these columns should be labeled "document_id" and "date", but the column labels can also be set using the `id.var` and `date.var` parameters.
#' Any other columns will automatically be included as document meta information in the output. 
#' 
#' The calculation of document similarity is performed using a vector space model approach. 
#' Inner-product based similarity measures are used, such as cosine similarity.
#' It is recommended to weight the DTM beforehand, for instance using Term frequency-inverse document frequency (tf.idf)
#' 
#' @param dtm A document-term matrix in the tm \link[tm]{DocumentTermMatrix} class. It is recommended to weight the DTM beforehand, for instance using \link[tm]{weightTfIdf}.
#' @param meta A data.frame where rows are documents and columns are document meta information. 
#' Should at least contain 2 columns: the document name/id and date. 
#' The name/id column should match the document names/ids of the edgelist, and its label is specified in the `id.var` argument. 
#' The date column should be intepretable with \link[base]{as.POSIXct}, and its label is specified in the `date.var` argument.            
#' @param id.var The label for the document name/id column in the `meta` data.frame. Default is "document_id"
#' @param date.var The label for the document date column in the `meta` data.frame . default is "date"
#' @param hour.window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. For example, c(-10, 36) will compare each document to all documents between the previous 10 and the next 36 hours.
#' @param measure the measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine" (cosine similarity), and the assymetrical measures "overlap_pct" (percentage of term scores in the document that also occur in the other document).
#' @param min.similarity a threshold for similarity. lower values are deleted. Set to 0.1 by default.
#' @param n.topsim An alternative or additional sort of threshold for similarity. Only keep the [n.topsim] highest similarity scores for x. Can return more than [n.topsim] similarity scores in the case of duplicate similarities.
#' @param only.from A vector with names/ids of documents (dtm rownames), or a logical vector that matches the rows of the dtm. Use to compare only these documents to other documents. 
#' @param only.to A vector with names/ids of documents (dtm rownames), or a logical vector that matches the rows of the dtm. Use to compare other documents to only these documents.
#' @param return.zeros If true, all comparison results are returned, including those with zero similarity (rarely usefull and problematic with large data)
#' @param only.complete.window if True, only compare articles (x) of which a full window of reference articles (y) is available. Thus, for the first and last [window.size] days, there will be no results for x.
#' 
#' @return A network/graph in the \link[igraph]{igraph} class
#' @export
#' 
#' @examples 
#' data(dtm)
#' data(meta)
#' 
#' dtm = tm::weightTfIdf(dtm)
#' g = newsflow.compare(dtm, meta, hour.window = c(0.1, 36))
#' 
#' vcount(g) # number of documents, or vertices
#' ecount(g) # number of document pairs, or edges
#' 
#' head(igraph::get.data.frame(g, 'vertices'))
#' head(igraph::get.data.frame(g, 'edges'))
newsflow.compare <- function(dtm, meta, id.var='document_id', date.var='date', hour.window=c(-24,24), measure='cosine', min.similarity=0, n.topsim=NULL, only.from=NULL, only.to=NULL, return.zeros=FALSE, only.complete.window=TRUE){
  confirm.dtm.meta(meta, id.var, date.var)
  meta = match.dtm.meta(dtm, meta, id.var)
  
  message('Indexing articles by date/time')
  if(is.null(only.from) & is.null(only.to)){
    dateids.x = dateids.y = getDateIds(meta[,date.var])
  } else{ 
    if(is.null(only.from)) only.from = rep(TRUE, nrow(dtm))
    if(is.null(only.to)) only.to = rep(TRUE, nrow(dtm))
    if(!class(only.from) == 'logical') only.from = rownames(dtm) %in% only.from
    if(!class(only.to) == 'logical') only.to = rownames(dtm) %in% only.to
    dateids.x = getDateIds(meta[,date.var], only.from)
    dateids.y = getDateIds(meta[,date.var], only.to)
  }
  dateindex = which(lapply(dateids.x, length) > 0)
  
  window = floor(hour.window[1]/24):ceiling(hour.window[2]/24)
  if(only.complete.window){
    if(window[1] < 0) dateindex = dateindex[dateindex > window[1]]
    if(rev(window)[1] > 0) dateindex = dateindex[dateindex <= length(dateids.x) - rev(window)[1]]  
  }
  
  message('Comparing documents')
  output = plyr::ldply(dateindex, function(i) ldply_documents.compare(i, dtm, dateids.x, dateids.y, window, measure, min.similarity, n.topsim, return.zeros), .progress='text')
  output = output[,!colnames(output) == '.id']
  
  message('Matching document meta')
  g = document.network(output, meta, id.var, date.var)
  
  delete.pairs = which(igraph::E(g)$hourdiff < hour.window[1] | igraph::E(g)$hourdiff > hour.window[2])
  g = igraph::delete.edges(g, delete.pairs)
  g
}

ldply_documents.compare <- function(i, dtm, dateids.x, dateids.y, window, measure, min.similarity, n.topsim, return.zeros){
  ## special function to be used in ldply in document.window.compare
  dtm.x_indices = unique(dateids.x[[i]])
  dtm.y_indices = unique(unlistWindow(dateids.y,i,window))
  if(length(dtm.y_indices) == 0) return(NULL)

  documents.compare(dtm[dtm.x_indices,], dtm[dtm.y_indices,], measure, min.similarity, n.topsim, return.zeros)
}


#' Delete duplicate (or similar) documents from a document term matrix 
#' 
#' Delete duplicate (or similar) documents from a document term matrix. 
#' Duplicates are defined by: having high content similarity, occuring within a given time distance and being published by the same source.
#' 
#' Note that this can also be used to delete "updates" of articles (e.g., on news sites, news agencies). 
#' This should be considered if the temporal order of publications is relevant for the analysis. 
#' 
#' @param dtm A document-term matrix in the tm \link[tm]{DocumentTermMatrix} class. It is recommended to weight the DTM beforehand, for instance using \link[tm]{weightTfIdf}.
#' @param meta A data.frame where rows are documents and columns are document meta information. 
#' Should contain 3 columns: the document name/id, date and source. 
#' The name/id column should match the document names/ids of the edgelist, and its label is specified in the `id.var` argument. 
#' The date column should be intepretable with \link[base]{as.POSIXct}, and its label is specified in the `date.var` argument.            
#' The source column is specified in the `date.var` argument.  
#' @param id.var The label for the document name/id column in the `meta` data.frame. Default is "document_id"
#' @param date.var The label for the document date column in the `meta` data.frame . default is "date"
#' @param source.var The label for the document date column in the `meta` data.frame . default is "source"
#' @param hour.window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. By default c(-24,24), which compares each document to all other documents within a 24 hour time distance.
#' @param measure the measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine" (cosine similarity), and the assymetrical measures "overlap_pct" (percentage of term scores in the document that also occur in the other document).
#' @param similarity a threshold for similarity. Documents of which similarity is equal or higher are deleted
#' @param keep A character indicating whether to keep the 'first' or 'last' published of duplicate documents.
#' @param tf.idf if TRUE, weight the dtm with tf.idf before comparing documents. The original (non-weighted) DTM is returned.
#' 
#' @return A dtm with the duplicate documents deleted
#' @export
#' 
#' @examples
#' data(dtm)
#' data(meta)
#' 
#' ## example with very low similarity threshold (normally not recommended!)
#' dtm2 = delete.duplicates(dtm, meta, similarity = 0.5, keep='first', tf.idf = TRUE)
delete.duplicates <- function(dtm, meta, id.var='document_id', date.var='date', source.var='source', hour.window=c(-24,24), measure='cosine', similarity=1, keep='first', tf.idf=FALSE){
  if(tf.idf) {
    g = newsflow.compare(tm::weightTfIdf(dtm), meta, measure=measure, min.similarity = similarity, hour.window=hour.window)
  } else {
    g = newsflow.compare(dtm, meta, measure=measure, min.similarity = similarity, hour.window=hour.window)
  }
  
  e = igraph::get.edges(g, igraph::E(g))
  d = igraph::get.data.frame(g, 'edges')  
  d$med.x = igraph::V(g)$source[e[,1]]
  d$med.y = igraph::V(g)$source[e[,2]]
  d = d[d$med.x == d$med.y,]
  
  duplicates = c()
  if(keep == 'first') {
    duplicates = c(duplicates, as.character(unique(d$to[d$hourdiff > 0])))
    duplicates = c(duplicates, as.character(unique(d$from[d$hourdiff < 0])))
  }
  if(keep == 'last') {
    duplicates = c(duplicates, as.character(unique(d$from[d$hourdiff > 0])))
    duplicates = c(duplicates, as.character(unique(d$to[d$hourdiff < 0])))
  }
  d = d[!d$from %in% duplicates & !d$to %in% duplicates,]
  
  ## if there are identical articles that occured simultaneously, delete randomly
  d = d[sample(1:nrow(d), nrow(d)),]
  d$fromi = match(d$from, unique(d$from, d$to))
  d$toi = match(d$to, unique(d$from, d$to))
  d = d[d$fromi < d$toi,]
  duplicates = unique(c(duplicates, as.character(d$from)))
  
  message('Deleting ', length(duplicates), ' duplicates')
 
  duplicates.med = meta$source[match(duplicates, rownames(dtm))]
  counts.med = table(duplicates.med)
  for(source in names(counts.med)){
      message('\t',source, ': ', counts.med[source])
  }
  
  dtm[!rownames(dtm) %in% duplicates,]
}

############

