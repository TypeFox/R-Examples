collapse <-
function(resp,group,minCell) {
    freq.table<-table(group,resp)
    if (ncol(freq.table)<2) stop("valid response categories must be >= 2")
    if (nrow(freq.table)<2) stop("valid group levels must be >= 2")
    original<-sort(unique(resp))
    suff.cell<-freq.table>=minCell
    suff.cat<-apply(suff.cell,2,prod)
    modified<-cumsum(suff.cat)
    modified[modified==0]<-1
    if (modified[1]==0) modified<-modified+1
    if (max(modified)<2) stop(paste("items must have at least two valid response categories with",minCell,"or more cases."))
    return (recode(resp,original,modified))
  }
