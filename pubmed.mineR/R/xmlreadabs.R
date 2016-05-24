xmlreadabs = function(file)
{test1 = readLines(file)
test1a = regexpr("PubmedArticle", test1, fixed=TRUE )
test1b = which(test1a != -1) 
test1c = lapply(seq(3,(length(test1b)-1),2), function(x){tempAA = test1b[x]; tempBB = test1b[(x+1)]; paste(test1[tempAA:tempBB], sep = "\n", collapse="")})
test1d = lapply(test1c, function(x){xmlTreeParse(x, asText=TRUE, useInternalNodes=TRUE)})
test1e = lapply(test1d, function(x){tempAA = unlist(getNodeSet(x, "//Abstract"));  tempBB = unlist(lapply(tempAA, function(x){xmlValue(x)})); if (length(tempBB) == 0) tempBB = "No Abstract Found"; return(tempBB[1]) }) 
test1f = lapply(test1d, function(x){tempAA = unlist(getNodeSet(x, "//ISOAbbreviation"));  tempBB = unlist(lapply(tempAA, function(x){xmlValue(x)})); if (length(tempBB) == 0 ) tempBB = "No Journal Found"; return(tempBB[1]) }) 
test1g = lapply(test1d, function(x){tempAA = unlist(getNodeSet(x, "//PMID"));  tempBB = unlist(lapply(tempAA, function(x){xmlValue(x)})); if (length(tempBB) == 0 ) tempBB = "No PMID Found"; return(tempBB[1])}) 
check = (length(test1e) == length(test1f)) &  (length(test1f) == length(test1g))

if (check) {resultabs = new("Abstracts",Journal = unlist(test1f), Abstract = unlist(test1e), PMID = as.numeric(unlist(test1g)));return(resultabs)}   else return("There is some problem in xml file. Please check") }
