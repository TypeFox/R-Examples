pubtator_function = function (x) 
{
    test = getURL(paste("http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/PubTator/abstract_ann.cgi?Disease=1&Gene=1&Chemical=1&Mutation=1&Species=1&pmid=", 
        x, sep = ""))
    testa = unlist(strsplit(test, "\n", fixed = T))
    table1 = NULL
    for (i in 4:length(testa)) {
	temps = unlist(strsplit(testa[i], "\t", fixed = T))
	if(length(temps) == 5) {temps = c(temps,"No Data")}
        table1 = rbind(table1, temps)
    }
    if (ncol(table1) == 6) {
        table2 = table1
        colnames(table2) = c("PMID", "Start", "End", "Term", 
            "TermType", "TermID")
        gene = NULL
        disease = NULL
        mutation = NULL
        chemical = NULL
        species = NULL
        for (i in 1:length(table2[, 5])) {
            if (table2[i, 5] == "Gene") 
                gene = c(gene, table2[i, 4])
            else if (table2[i, 5] == "Disease") 
                disease = c(disease, table2[i, 4])
            else if (table2[i, 5] == "Mutation") 
                mutation = c(mutation, table2[i, 4])
            else if (table2[i, 5] == "Chemical") 
                chemical = c(chemical, table2[i, 4])
            else if (table2[i, 5] == "Species") 
                species = c(species, table2[i, 4])
        }
        gene = union(gene, gene)
        disease = union(disease, disease)
        mutation = union(mutation, mutation)
        chemical = union(chemical, chemical)
        species = union(species, species)
        return(list(Genes = gene, Diseases = disease, Mutations = mutation, 
            Chemicals = chemical, Species = species, PMID = x))
    }
    else return(" No Data ")
}
