metalist.to.matrix<-function(list,genenames=NULL)
{
n.study<-which(names(list)=="TestStatistic")
xx<-matrix(0,nrow=length(list$TestStatistic), ncol=n.study)
for (i in 1:(n.study-1))
{xx[list[[i]],i]<-1}
xx[,n.study]<-list[[n.study]]
colnames(xx)<-names(list)[1:n.study]
if (is.null(genenames)) rownames(xx)<-list[["gene.names"]] else rownames(xx)<-genenames
return(xx)
}

