access_import <-
function(file,table_names,ntab=length(table_names),SQL_query=rep(F,times=ntab),where_sql=c(),out.format="d-m-yy",uid="",pwd="",...)
{
  mycon<-RODBC::odbcConnectAccess(file,uid=uid,pwd=pwd,...)
  for(i in 1:length(where_sql)){SQL_query[where_sql[i]]<-T}
  if(ntab==1){
    if(!SQL_query){x<-sqlFetch(mycon,table_names,...);close(mycon)
    for(i in 1:length(x[1,])){
      if(class(x[1,i])[1]==c("POSIXct")){x[,i]<-chron(as.character(x[,i]),format=("y-m-d"),out.format=out.format,...)}}
    return(x)}
    else{x<-sqlQuery(mycon,table_names,SQL_query,...);close(mycon)
    for(i in 1:length(x[1,])){
      if(class(x[1,i])[1]==c("POSIXct")){x[,i]<-chron(as.character(x[,i]),format=("y-m-d"),out.format=out.format,...)}}
    return(x)}
  }
  x<-vector("list",ntab)
  for(i in 1:ntab){
    if(!SQL_query[i]){x[[i]]<-sqlFetch(mycon,table_names[i],...)
    for(j in 1:length(x[[i]][1,])){
      if(class(x[[i]][1,j])[1]==c("POSIXct")){x[[i]][,j]<-chron(as.character(x[[i]][,j]),format=("y-m-d"),out.format=out.format,...)}}
    }
    else{x[[i]]<-sqlQuery(mycon,table_names[i],...)
    for(j in 1:length(x[[i]][1,])){
      if(class(x[[i]][1,j])[1]==c("POSIXct")){x[[i]][,j]<-chron(as.character(x[[i]][,j]),format=("y-m-d"),out.format=out.format,...)}}
    }
  }
  close(mycon)
  names(x)<-table_names
  return(x)
}
