		#get a list of coeffient for each treatment
getTrtCoef <-
function(design.df, trtCols){

if(length(trtCols) == 1){
return(length(design.df[,trtCols])/nlevels(design.df[,trtCols]))
}else if(any(grepl(":", trtCols))){ 

trtColsList = lapply(strsplit(trtCols, "\\:"), function(x) design.df[,x])
repList = sapply(trtColsList, function(y) if(is.factor(y)) {mean(table(y))} 
else{ mean(table(apply(y, 1, function(x) paste(x, collapse ="."))))} )

return(repList)      
}else{
repList = apply(design.df[,trtCols], 2, function(x)mean(table(x)))
return(repList)      
}                                                                 
}
