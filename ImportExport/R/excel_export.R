excel_export <-
function(x,file,table_names=as.character(1:length(x)),row.names=F,...) {
  
  if(class(x)=="data.frame"){xlsx::write.xlsx(x,file,sheetName=table_names[1],row.names=row.names,...)}
  if(class(x)=="list"){ for(i in 1:length(x)){ xlsx::write.xlsx(x[[i]],file,sheetName=table_names[i],append=(i!=1),row.names=row.names,...)}}
}
