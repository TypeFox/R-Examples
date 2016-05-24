if(getRversion() >= "2.15.0") utils::globalVariables("common_words_new");
word_atomizations = function(m){tempzz = unlist(lapply(m@Abstract, function(x){tempa = strsplit(x, ".  ",fixed=T);tempa1 = which( nchar(tempa[[1]]) == max(nchar(tempa[[1]])));
tempb = unlist(strsplit(tempa[[1]][tempa1], ".",fixed = T));
tempc = unlist(strsplit(tempb, ",",fixed = T));
tempd = unlist(  strsplit(tempc, ":",fixed = T));
tempe = unlist( strsplit(tempd, ";",fixed = T));
tempe1 = unlist( strsplit(tempe, "'",fixed = T));
tempf = unlist(  strsplit(tempe1, " ",fixed = T));
tempg = tolower(tempf);
temph = sort(tempg);
return(tempg)}));
tempi = as.data.frame(table(tempzz));
tempj = unlist(lapply(common_words_new, function(x){tempoo = which(as.character(tempi[,1]) == x   ); if (length(tempoo) != 0) return(tempoo)}));
tempk = tempi[-tempj,];
tempk2=tempk[order(as.numeric(tempk$Freq), decreasing = T),]
colnames(tempk2)=c("words", "Freq")
write.table(tempk2, file = "word_table.txt", sep = "\t", row.names = F); return(tempk2)}
