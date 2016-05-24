printabs <- function(object){ cat(paste("Number of Abstracts",length(object@PMID)), fill = T)  ; 
tempa <- SentenceToken(object@Abstract[1]);
tempb <- length(object@PMID); 
tempc <-  SentenceToken(object@Abstract[tempb]);   
cat("Starts with ", fill = T); 
cat(tempa[1:3], fill = T);  
cat("  ", fill = T); 
cat("                          ", fill = T); 
cat("Ends with ",fill = T);
cat(tempc[1:3], fill = T)   }
