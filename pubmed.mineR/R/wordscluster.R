wordscluster = function(lower,upper)
{ 

words = readLines("word_table.txt")

checkword = strsplit(words,"\"", fixed=T)

checkword2 = unlist(lapply(checkword, function(x){return(x[2])}))


  if (missing(lower)) lower = 5;
  if (missing(upper)) upper = 30;

wordindex1 = which(nchar(checkword2) >= lower )
subcheckword2_1 = checkword2[wordindex1]

wordindex2 = which(nchar(subcheckword2_1) <= upper )
subcheckword2_2 = subcheckword2_1[wordindex2]

filteredsubcheckword2_2 = unlist(lapply(subcheckword2_2, function(x){tempwww = unlist(strsplit(x,"",fixed=T)); tempxxx =  isTRUE(all.equal(tempwww %in% c(letters,"-"), rep(TRUE, length(tempwww)))) ; if (tempxxx) return(x)}))
range(nchar(filteredsubcheckword2_2))
ncharfilteredsubcheckword2_2 = nchar(filteredsubcheckword2_2)
ncharfilteredsubcheckword2_2sorted = sort(ncharfilteredsubcheckword2_2, index.return = T)
tempaaa = ncharfilteredsubcheckword2_2sorted$ix
filteredsubcheckword2_2sorted = filteredsubcheckword2_2[tempaaa]

newwordscluster = NULL;  tempC = filteredsubcheckword2_2sorted; repeat  {print(paste("Number of words remaining", length(tempC), "wordsclustersize", length(newwordscluster), sep = " "));  if ( length(tempC) > 200) xa = 200 else xa = length(tempC);    tempB = unlist(lapply(tempC[1:xa], function(x){tempA = agrep(x, tempC, max.distance = 0, value = F); return(tempA)}));   if (length(tempB) != 0 )  {newwordscluster =  c(newwordscluster,   lapply(tempC[1:xa], function(x){tempA = agrep(x, tempC, max.distance = 0, value = T); return(tempA)})); tempB = union(tempB,tempB); tempC = tempC[-tempB]} else break}

forsorting = sort(unlist(lapply(newwordscluster, function(x){return(x[1])}  )), index.return=T)
newwordsclustersorted  = lapply(forsorting$ix, function(x){return(unlist(newwordscluster[[x]]))}) 
repword = unlist(lapply(newwordsclustersorted, function(x){return(x[1])}  ))
newwordsclustersize = unlist(lapply(newwordsclustersorted, function(x){return(length(x))}))
resulttable = cbind(newwordsclustersize, repword)
colnames(resulttable) = c("Cluster Size", "Representative Word")
write.table(resulttable, file = "resulttable.txt", quote = F, sep = "\t"); return(newwordsclustersorted)}