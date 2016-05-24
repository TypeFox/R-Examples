"CVratio" <-
function(NTOTAL, NESSENTIAL)
 {
 n <- NTOTAL
 ne <- NESSENTIAL
 cvr <- (ne - n/2)/(n/2)
 return(cvr)
 }

