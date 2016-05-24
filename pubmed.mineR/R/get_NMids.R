get_NMids = function(x){
test1 = getURL(paste("http://www.ncbi.nlm.nih.gov/gene/?term=",x,sep=""))
test1a = htmlTreeParse(test1,useInternalNodes = TRUE)
test1b = lapply(getNodeSet(test1a,"//a"),function(x){xmlValue(x)})
resA = NULL; for (i in 1:length(test1b)){temp = test1b[[i]];temp1 = regexpr("NM_",temp);if (temp1!= -1)  resA = c(resA,test1b[[i]])}; return(resA)}
