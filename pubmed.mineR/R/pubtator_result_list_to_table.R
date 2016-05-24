pubtator_result_list_to_table=function (x) 
{
    write(c("Genes", "Diseases", "Mutations", "Chemicals", "Species", 
        "PMID"), file = "result.txt", append = T, sep = "\t", 
        ncolumns = 6)
    for (i in 1:length(x)) {
        if (x[[i]][1] != " No Data ") 
            write.table(cbind(paste(x[[i]]$Genes, collapse = ","), 
                paste(x[[i]]$Diseases, collapse = ","), paste(x[[i]]$Mutations, 
                  collapse = ","), paste(x[[i]]$Chemicals, collapse = ","), 
                paste(x[[i]]$Species, collapse = ","), x[[i]]$PMID), file = "result.txt", 
                sep = "\t", quote = F, row.names = F, col.names = F, 
                append = T)
    }
}

