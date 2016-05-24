Give_Sentences_PMC = function(PMCID,term){
  check = getURL(paste("http://www.ncbi.nlm.nih.gov/pmc/articles/PMC", PMCID, "/", sep=""))
                 checkxml = htmlTreeParse(check,useInternalNodes = T)
                 check3 = lapply(getNodeSet(checkxml,"//p"),function(x){xmlValue(x)})
                 test = NULL; for (i in 1:length(check3)) {tempA = SentenceToken(check3[[i]]);tempB = regexpr(term,tempA); tempC = which(attr(tempB,"match.length") == nchar(term)); if ( length(tempC) != 0 ) tempD = tempA[tempC] else tempD = NULL; test = c(test,tempD)}; return(test)}
