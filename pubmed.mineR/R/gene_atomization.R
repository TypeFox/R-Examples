if(getRversion() >= "2.15.0") utils::globalVariables(c("common_words_new", "HGNCdata"));
gene_atomization <-function(m){tempzz = unlist(lapply(m@Abstract, function(x){tempa = strsplit(x, ".  ",fixed=T);tempa1 = which( nchar(tempa[[1]]) == max(nchar(tempa[[1]])));
tempb = unlist(strsplit(tempa[[1]][tempa1], ".",fixed = T));
tempc = unlist(strsplit(tempb, ",",fixed = T));
tempd = unlist(  strsplit(tempc, ":",fixed = T));
tempe = unlist( strsplit(tempd, ";",fixed = T));
tempe1 = unlist( strsplit(tempe, "'",fixed = T));
tempf = unlist(  strsplit(tempe1, " ",fixed = T));
return(tempf)}));
tempi = as.data.frame(table(tempzz));
tempj = unlist(lapply(common_words_new, function(x){tempoo = which(as.character(tempi[,1]) == x   ); if (length(tempoo) != 0) return(tempoo)}));
tempk = tempi[-tempj,];
templ = as.character(HGNCdata$Approved.Symbol);
tempm = unlist(lapply(templ,function(x){return(which(x == as.character(tempk$tempzz) ))}));
tempn = tempk[tempm,]
tempn2=tempn[order(as.numeric(tempn$Freq), decreasing = T),]
tempo = unlist(lapply( as.character(tempn2$tempzz), function(x){return(which(x == templ))}));
Genes = as.character(HGNCdata$Approved.Name[tempo]);
data_table= cbind(as.character(tempn2$tempzz), Genes, tempn2$Freq);
colnames(data_table)=c("Gene_symbol", "   Genes", "Freq");
write.table(data_table, file = "table.txt", sep = "\t", row.names = F); return(data_table) }
