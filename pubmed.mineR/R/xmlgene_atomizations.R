if(getRversion() >= "2.15.0") utils::globalVariables(c("common_words_new", "HGNCdata"));
xmlgene_atomizations = function (m) 
{
    tempb = unlist(lapply(m@Abstract, function(x) {tempa = strsplit(x, ".", fixed = T); return(tempa)}))
   tempc = unlist(strsplit(tempb, ",", fixed = T))
   tempd = unlist(strsplit(tempc, ":", fixed = T))
   tempe = unlist(strsplit(tempd, ";", fixed = T))
   tempe1 = unlist(strsplit(tempe, "'", fixed = T))
   tempf = unlist(strsplit(tempe1, " ", fixed = T)) 
    tempi = as.data.frame(table(tempf))
    tempj = unlist(lapply(common_words_new, function(x) {
        tempoo = which(as.character(tempi[, 1]) == x)
        if (length(tempoo) != 0) 
            return(tempoo)
    }))
    tempk = tempi[-tempj, ]
    templ = as.character(HGNCdata$Approved.Symbol)
    tempm = unlist(lapply(templ, function(x) {
        return(which(x == as.character(tempk$tempf)))
    }))
    tempn = tempk[tempm, ]
    tempn2 = tempn[order(as.numeric(tempn$Freq), decreasing = T), 
        ]
    tempo = unlist(lapply(as.character(tempn2$tempf), function(x) {
        return(which(x == templ))
    }))
    Genes = as.character(HGNCdata$Approved.Name[tempo])
    data_table = cbind(as.character(tempn2$tempf), Genes, tempn2$Freq)
    colnames(data_table) = c("Gene_symbol", "   Genes", "Freq")
    write.table(data_table, file = "table.txt", sep = "\t", row.names = F)
    return(data_table)
}
