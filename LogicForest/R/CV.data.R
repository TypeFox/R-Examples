CV.data <-
function(resp, Xs, kfold)
{
 data<-cbind(Xs, resp)
 n<-nrow(Xs)
 sz<-floor(n/kfold)
 sz2<-sz+1
 rmdr<-n-kfold*sz
 sz.vec<-rep(c(sz2, sz), c(rmdr, kfold-rmdr))
 sets<-vector("list", kfold)
 for (i in 1:kfold)
   {
   len<-sz.vec[i]
   len.sum<-sum(sz.vec[1:i])
   if(i==1) {
      sub<-c(1:len)
      testdata<-data[sub,]
      traindata<-data[-sub,]
      }
   if(i>1) {
      sub<-c(((len.sum-len)+1):len.sum)
      testdata<-data[sub,]
      traindata<-data[-sub,]
      }
   sets[[i]]<-list(traindata=traindata, testdata=testdata, loc.ids=sub)
   }
 sets
}
