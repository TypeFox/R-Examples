filter.site<-function(data,site)
{
  if(!is.data.frame(data) || !is.character(site) || !grepl("^[a-zA-Z]+$",gsub(" ", "", site)))
     stop("invalid input parameter(s) specification: check data/site")
  
  data(vmdbsite,envir=environment())
  vmdbsite<-get("vmdbsite",envir=environment())
  site.code<-vmdbsite[substr(as.character(vmdbsite$Name),start=1,stop=nchar(site))==toupper(site),1]
  data[data$sitecode==site.code,]
}
