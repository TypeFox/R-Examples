ARorder <-
function(x,pmax,type)
{
x<-as.matrix(x)
if (type=="const")
M <- AR.order(x,pmax)
if (type=="const+trend")
M <- ART.order(x,pmax)
rownames(M$Criteria) <- paste("",1:pmax,sep="") 
rownames(M$ARorder) <- "p*"; colnames(M$ARorder) <- c("aic","bic","hq")
return(list(ARorder=M$ARorder,Criteria=M$Criteria))
}
