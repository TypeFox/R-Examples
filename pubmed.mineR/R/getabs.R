setGeneric("getabs", function(object,x,y) standardGeneric("getabs"));
setMethod("getabs","Abstracts",function(object,x,y){
tempabs <- object@Abstract; tempuids <- object@PMID; tempjournal <- object@Journal; 
xj<- regexpr(x,tempabs,ignore.case = y); xxj<- 0;xxj<- c(xxj,which(xj != -1)); 
newtempuids <- NULL; newtempabs <- NULL; newtempjournal <- NULL; 
for(m in seq(along = xxj)) {if(xxj[m] > 0 ) {tempm <- xxj[m];newtempuids<- c(newtempuids,tempuids[tempm]); newtempabs <- c(newtempabs,tempabs[tempm]); newtempjournal <- c(newtempjournal,tempjournal[tempm]) }}; 
 if ( length(xxj) == 1 ) {print(paste("No abstracts found for the term",x,sep =" ")); write(paste("No abstracts found for the term",x,sep =" "), file = "dataout.txt", append = T); result <- new("Abstracts", Journal = "NONE", Abstract = "NONE", PMID = 0); return(result)} else {result <- new("Abstracts", Journal = newtempjournal, Abstract = newtempabs, PMID = newtempuids); print(paste(length(result@PMID),"abstracts",x,sep = " "));   write(paste(length(result@PMID),"abstracts",x,sep = " "), file = "dataout.txt", append = T);return(result)}  })
