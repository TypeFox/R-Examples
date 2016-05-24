`extractPval` <-
function(x)
 {
   
   models<-attr(x,"models")
   if(length(models)==6)
    models<-c(1:5) 
   quantitative<-attr(x,"quantitative")
   pos<-ifelse(quantitative,7,8)

   ans<-t(data.frame(lapply(1:length(x),extractPval.i,x=x,pos=pos,models=models)))

   ans<-data.frame(ans)
   for (i in 2:ncol(ans))
      ans[,i]<-as.numeric(as.character(ans[,i]))

   dimnames(ans)[[1]]<-attr(x,"label.SNPs")
   dimnames(ans)[[2]]<-c("comments",c("codominant","dominant","recessive","overdominant","log-additive")[models])
   ans

 }

