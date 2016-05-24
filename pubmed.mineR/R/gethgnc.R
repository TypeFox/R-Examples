setGeneric("gethgnc", function(object,x) standardGeneric("gethgnc"));
setMethod("gethgnc","HGNC",function(object,x){ check1 <- match(x, object@ApprovedSymbol, nomatch = 0); 
check2 <- regexpr(x, object@ApprovedName, ignore.case=TRUE); check3 <-  regexpr(x, object@Aliases, ignore.case=TRUE); 
xxa <- NULL; if (check1 !=0) {xxa <- c(object@ApprovedName[check1],object@Aliases[check1])}; 
for (j in 1:19363 ) { if (check2[j] != -1) xxa <- c(xxa,object@ApprovedSymbol[j],object@Aliases[j])  else  if (check3[j] != -1) xxa <- c(xxa,object@ApprovedSymbol[j],object@ApprovedName[j])};  return(xxa) })
