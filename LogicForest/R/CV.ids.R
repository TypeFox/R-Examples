CV.ids <-
function(n, kfold)
{
 sz<-floor(n/kfold)
 sz2<-sz+1
 rmdr<-n-kfold*sz
 sz.vec<-rep(c(sz2, sz), c(rmdr, kfold-rmdr))
 loc.ids<-vector("list", kfold)
 for (i in 1:kfold)
   {
   len<-sz.vec[i]
   len.sum<-sum(sz.vec[1:i])
   if(i==1) {sub<-c(1:len)}
   if(i>1) {sub<-c(((len.sum-len)+1):len.sum)}
   loc.ids[[i]]<-sub
   }
 loc.ids
}
