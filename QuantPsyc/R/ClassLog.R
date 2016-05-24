"ClassLog" <- function(MOD, resp, cut=.5)
{
 classtab <- prop.table (table (predict(MOD, type='response') > cut, resp), 2)
 rawtab <- table (predict(MOD, type='response') > cut, resp)
 overall <- sum(diag(prop.table (table (predict(MOD, type='response') > cut, resp))))
 mcFadden <- 1 - MOD$deviance/MOD$null.deviance 
 return(
 	list(
 		rawtab=rawtab, 
 		classtab=classtab, 
 		overall=overall, 
 		mcFadden=mcFadden)
 		)
 }


