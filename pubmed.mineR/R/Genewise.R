setGeneric("Genewise", function(object,gene) standardGeneric("Genewise"));
setMethod("Genewise", "Abstracts", function(object, gene)
{
tempj=object@Abstract;
tempk=regexpr(as.character(gene), tempj, ignore.case=FALSE);
templ=which(tempk != -1);
if(length(templ)==0)
{print(paste("No abstracts found for the term",gene,sep =" ")); 
write(paste("No abstracts found for the term",gene,sep =" "), file = "dataout.txt", append = T)}
else {print(paste(length(templ),"abstracts",gene, sep = " "));
write(paste(length(templ),"abstracts",gene, sep = " "), file = "dataout.txt", append = T)} })
