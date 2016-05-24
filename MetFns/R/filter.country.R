filter.country<-function(data,country)
{
  if(!is.data.frame(data) || !is.character(country) || !grepl("^[a-zA-Z]+$",gsub(" ", "", country)))
     stop("invalid input parameter(s) specification: check data/country")
 
  data(vmdbsite,envir=environment())
  vmdbsite<-get("vmdbsite",envir=environment())
  country.codes<-vmdbsite[as.character(vmdbsite$Country)==toupper(country),1]
  data[which(data$sitecode%in%country.codes),]
}
