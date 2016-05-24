mloglike_offsetpois <-
function(data1, para)
{
   N=length(data1)   
   data11<-subset(data1, data1>=para[4])
   data12<-subset(data1, data1<para[4])
   temp=log(para[1])+dpois(data11-para[4], para[2], log=TRUE)
   temp1=log(1+exp(log(1-para[1])+dpois(data11, para[3], log=TRUE)-log(para[1])-dpois(data11-para[4],para[2], log=TRUE)))
   temp2=dpois(data12, para[3], log=TRUE)
   loglike=sum(temp)+sum(temp1)+sum(temp2)
   return(loglike)
}
