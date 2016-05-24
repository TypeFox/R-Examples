spss_import <-
function(file, allow="_",out.format="d-m-yy",use.value.labels=F,...){
  varlist<-data.frame(spss_varlist(file))
  vardate<-subset(varlist,substring(varlist$printfmt,1,4)=="DATE")$longname
  x<-try((Hmisc::spss.get(file,datevars=vardate,allow=allow,use.value.labels=use.value.labels,...)),silent=T)
  if(class(x)=="data.frame"){
  for (i in vardate) { descr<-attr(x[,i],"label")
    x[,i]<-chron(as.character(x[,i]),format="y-m-d",out.format=out.format)
    attr(x[,i],"label")<-descr}}
  else{ x<-read_sav(file);cat("readed with read_sav")}
  
  x
}
