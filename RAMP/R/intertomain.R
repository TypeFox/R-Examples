intertomain <-
function(interNameList, p){
###convert an interaction name list into the corresponding main effect for enforcing strong hierachy.
###input: interNameList: interaction name list
###input: p: dimensionality
mainInd = rep(0,p)
for(i in 1:length(interNameList)){
	interName = interNameList[[i]]
	pair = as.numeric(strsplit(interName, "X")[[1]][2:3])
	mainInd[pair[1]] = 1
	mainInd[pair[2]] = 1
}
return(which(mainInd==1))
}
