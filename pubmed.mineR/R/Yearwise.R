setGeneric("Yearwise", function(object,year) standardGeneric("Yearwise"));
setMethod("Yearwise", "Abstracts", function(object, year)
{
tempj=object@Journal;
tempk=regexpr(as.character(year), tempj, ignore.case=FALSE);
templ=which(tempk != -1);
if(length(templ)==0)
{print(paste("No abstracts found for the term",year,sep =" ")); 
write(paste("No abstracts found for the term",year,sep =" "), file = "dataout.txt", append = T)}
else {print(paste(length(templ),"abstracts",year, sep = " "));
write(paste(length(templ),"abstracts",year, sep = " "), file = "dataout.txt", append = T)} })
