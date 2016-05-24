"savecluster" <-
structure(function(v,filename,cont=FALSE){
if(length(grep(patt="w32",x=version["os"]))){
	eol<-"\n"
}else{eol<-"\r\n"}
cat(paste(v,collapse=eol),file = filename,append=cont)
}
, comment = "Save cluster to file that can be read by Pajek")
