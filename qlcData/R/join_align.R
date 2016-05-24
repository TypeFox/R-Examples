# =============================================
# join multialignments into one large dataframe
# =============================================

# duplicated are simply discarded
# names are aligned between different multialignments

join.align <- function(alignments) {
  
  names(alignments) <- NULL
  
  # which lects occur	
  lects <- sapply(alignments, function(x){x$doculects})
  unique_lects <- sort(unique(unlist(lects)))
  
  # reorder tables to join alignment into one table
  # and remove duplicated
  reorder <- function(align, l = unique_lects){
    a <- align$align
    d <- align$doculects
    dupl <- duplicated(d)
    a <- a[!dupl,]
    d <- d[!dupl]
    new_order <- as.numeric(factor(l, levels = d))
    return(a[new_order,])
  }
  ordered <- sapply(alignments, reorder)
  
  # join all alignments
  join <- do.call(cbind,ordered)
  rownames(join) <- unique_lects
  
  return(join)
  
}