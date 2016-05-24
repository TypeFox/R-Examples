filter.imocode<-function(data, imocode)
{
   if(!is.data.frame(data) || !is.character(imocode) ||  nchar(imocode)!=5 || !grepl("^[A-Z]+$",imocode))
      stop("invalid input parameter(s) specification: check data/imocode")
  
   data[data$IMOcode==imocode,]
}
