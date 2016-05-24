mloglike_offsetNB <-
function(data1, para)
{
   N=length(data1)   
   data11<-subset(data1, data1>=para[6])
   data12<-subset(data1, data1<para[6])
   temp=log(para[1])+dnbinom(data11-para[6], para[3],,para[2], log=TRUE)
   temp1=log(1+exp(log(1-para[1])+dnbinom(data11, para[5],,para[4], log=TRUE)-log(para[1])-dnbinom(data11-para[6], para[3],,para[2], log=TRUE)))
   temp2=dnbinom(data12, para[5], ,para[4], log=TRUE)
   loglike=sum(temp)+sum(temp1)+sum(temp2)
   return(loglike)
}
