statMissing<-
function (X, tri = 0, effectif=0, ...) 
{
    noms <- as.character(whatis(X)$variable.name)
    n<-nrow(X)
    
    if (tri == 0) {

if (effectif == 0){
manque <- 100*whatis(X)$missing/n
names(manque) <- paste(noms," (%)",sep="")
return(round(manque))}
	else{
manque <- whatis(X)$missing
names(manque) <- noms
return(manque)}
		
    		}
    
else{
if (effectif == 0){
manque <- 100*whatis(X)$missing/n
names(manque) <- paste(noms," (%)",sep="")
return(round(sort(manque)))}
else{
manque <- whatis(X)$missing
names(manque) <- noms
return(sort(manque))}
}
    
}