		#get the incidence matrix N
getIncidenceMatrix <-
function(design.df, trtCols){

if(length(trtCols) == 1){
incident = design.df[,trtCols]
nLevels = levels(design.df[,trtCols])    
}else if(any(grepl(":", trtCols))){ 
uniqueTrtCols = unique(unlist(strsplit(trtCols, "\\:")))  

incident = as.factor(apply(design.df[,uniqueTrtCols], 1, function(x) paste(x, collapse =".")))
nLevels = sort(levels(interaction(design.df[,uniqueTrtCols])))

}else{
incident = as.factor(apply(design.df[,trtCols], 1, function(x) paste(x, collapse =".")))
nLevels = sort(levels(interaction(design.df[,trtCols])))
}       

N <- matrix(0, nrow=nrow(design.df), ncol=length(nLevels))
N[cbind(1:nrow(design.df), match(incident, nLevels))] <- 1

return(N)
}
