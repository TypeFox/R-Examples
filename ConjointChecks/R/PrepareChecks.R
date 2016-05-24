PrepareChecks<-function(resp,ss.lower=10) {
  if (any(is.na(resp))) stop("Checks will only work with complete data. Suggestion: remove respondents with missing responses.")
  if (ss.lower==1) {
    message("ss.lower must be greater than 1, setting to 2.")
    ss.lower<-2
  }
  ncol(resp)->n.items
  #reorder things
  colSums(resp)->cs
  resp[,order(cs)]->resp
  #group by sum scores
  rowSums(resp)->rs
  n<-N<-list()
  table(rs)->tab
  as.numeric(names(tab))[tab>=ss.lower]->lev
  for (s in lev) {
    resp[rs==s,]->tmp
    rep(nrow(tmp),n.items)->N[[as.character(s)]]
    colSums(tmp)->n[[as.character(s)]]
  }
  do.call("rbind",N)->N
  do.call("rbind",n)->n
  list(N=N,n=n)
}
