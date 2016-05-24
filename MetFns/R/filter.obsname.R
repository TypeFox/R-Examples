filter.obsname<-function(data,name,fname)
{
   if(!is.data.frame(data) || !is.character(c(name,fname)) || !grepl("^[a-zA-Z]+$",name) || !grepl("^[a-zA-Z]+$",fname))
      stop("invalid input parameter(s) specification: check data/name/fname")
   
   data(vmdbpers,envir=environment())
   vmdbpers<-get("vmdbpers",envir=environment()) 
   vmdbpers <- data.frame(lapply(vmdbpers, as.character), stringsAsFactors=FALSE)
   imocode<-vmdbpers[vmdbpers$Name==toupper(name) & vmdbpers$Firstname==toupper(fname),1]
   filter.imocode(data,imocode)
}
