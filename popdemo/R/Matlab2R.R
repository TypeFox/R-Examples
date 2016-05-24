Matlab2R <-
function(M){
rows<-strsplit(M,";")[[1]]
rows<-strsplit(rows," ")
order<-length(rows)
for(i in 1:order){
    rows[[i]]<-paste(noquote(rows[[i]][!rows[[i]]==""]),collapse=",")
}
elements<-noquote(paste(rows,collapse=","))
elements<-gsub("\\[,","",elements)
elements<-gsub("\\[","",elements)
elements<-gsub("\\,]","",elements)
elements<-gsub("\\]","",elements)
Rmat<-character(0)
command<-paste("matrix(c(",elements,"),nrow=order,byrow=TRUE)")
Rmat<-eval(parse(text=noquote(command)))
return(Rmat)
}

