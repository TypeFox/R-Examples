read.write <-
function(rbind.result,link){
suppressWarnings(tryCatch(res<-read.table(link,header=T), error=function(e) write.table(matrix(c("ID","Question.number","type","Condition.Likert","Response"),nrow=1),link,quote=FALSE,row.names=FALSE,col.names=FALSE)))
res<-read.table(link,header=T)
res.new<-rbind(res,rbind.result)
write.table(res.new,paste(link),quote=FALSE,row.names=FALSE,col.names=TRUE)
}
