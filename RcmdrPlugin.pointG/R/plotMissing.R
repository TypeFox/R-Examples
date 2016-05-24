plotMissing<-
function (X, tri = 0, effectif=0, ...) 
{
    noms <- as.character(whatis(X)$variable.name)
    n<-nrow(X)
    
    if (tri == 0) {

if (effectif == 0){
manque <- 100*whatis(X)$missing/n
names(manque) <- noms
maximum <- max(manque)
        dotchart(manque, xlab = paste("Donn", "\U00E9", "es manquantes (%)", 
            sep = ""), xlim = c(0, maximum), ...)}

	else{
manque <- whatis(X)$missing
names(manque) <- noms
maximum <- max(manque)
dotchart(manque, xlab = paste("Donn", "\U00E9", "es manquantes", 
            sep = ""), xlim = c(0, maximum), ...)}
		
    		}
    
else{
if (effectif == 0){
manque <- 100*whatis(X)$missing/n
names(manque) <- noms
maximum <- max(manque)
        dotchart(sort(manque), xlab = paste("Donn", "\U00E9", "es manquantes (%)", 
            sep = ""), xlim = c(0, maximum), ...)}
else{
manque <- whatis(X)$missing
names(manque) <- noms
maximum <- max(manque)
dotchart(sort(manque), xlab = paste("Donn", "\U00E9", "es manquantes", 
            sep = ""), xlim = c(0, maximum), ...)}
}
    
}