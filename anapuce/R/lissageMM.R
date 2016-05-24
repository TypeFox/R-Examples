lissageMM <-
function(x,indice,mm) {
 
 tmp <- apply(as.matrix(c(1:5)),1,FUN=function(x,indice) unique(max(match(x,indice),0,na.rm=TRUE)),indice=indice)
 tmplist <- list()
 for (i in c(1:4)) {
 tmplist[[i]] <- fonctionMobile(x,mean,mm[i])
 }
 y<-x
 for (i in 1:length(x)){
    if ((i>=tmp[2])&&(i<tmp[3])) y[i]<-tmplist[[1]][i] 
    if ((i>=tmp[3])&&(i<tmp[4])) y[i]<-tmplist[[2]][i] 
    if ((i>=tmp[4])&&(i<tmp[5])) y[i]<-tmplist[[3]][i] 
    if (i>=tmp[5]) y[i]<-tmplist[[4]][i]
 }

return(y)
}

