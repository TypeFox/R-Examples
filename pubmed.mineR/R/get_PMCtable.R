get_PMCtable = function(url){
test = getURL(url)
test1 = htmlTreeParse(test,useInternalNodes = T)
test2 = lapply(getNodeSet(test1,"//tr"),function(x){(x)})
table=NULL;for (i in 3:length(test2)){table=rbind(table,getChildrenStrings(test2[[i]]))  }
colnames(table)= as.character(getChildrenStrings(test2[[2]]))
return(table)}
