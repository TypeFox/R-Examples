list.PIs <-
function(match.list, allPIs, preds)
{
  PIlist<-unique(match.list)
  PIsetlist<-vector("list", length(PIlist))
  allPI.list<-vector("list", length(allPIs))
  lgth<-c()
  all.lgth<-c()
  for (i in 1:length(PIlist))
    {
    PIstr<-unlist(strsplit(PIlist[[i]], " "))
    PIsetlist[[i]]<-setdiff(PIstr, "&")
    lg<-length(PIsetlist[[i]])
    lgth<-append(lgth, lg)
    }
  for (j in 1:length(allPIs))
    {  
    PIst<-unlist(strsplit(allPIs[[j]], " "))
    allPI.list[[j]]<-setdiff(PIst, "&")
    Alg<-length(allPI.list[[j]])
    all.lgth<-append(all.lgth, Alg)
    }
  lengths<-sort(setdiff(lgth, 0))
  all.lengths<-sort(setdiff(all.lgth, 0))
  ans<-list(PIsetlist=PIsetlist, allPI.list=allPI.list, lengths=lengths, all.lengths=all.lengths)
}
