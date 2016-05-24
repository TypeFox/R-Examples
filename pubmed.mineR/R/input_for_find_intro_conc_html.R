input_for_find_intro_conc_html = function (y, all) 
{
    check0 = lapply(y@PMID, function(b) {
        url = paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=", 
            "pubmed", "&id=", b, sep = "")
        epost = xmlTreeParse(getURL(url), useInternalNodes = T)
        webenv = xmlValue(getNodeSet(epost, "//WebEnv")[[1]])
        key = xmlValue(getNodeSet(epost, "//QueryKey")[[1]])
        url1 = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        query = "db=pubmed&retmode=xml&rettype=abstracts"
        efetch = xmlTreeParse(getURL(paste(url1, query, "&WebEnv=", 
            webenv, "&query_key=", key, sep = "")), useInternalNodes = T)
        abs = unlist(lapply(getNodeSet(efetch, "//Abstract"), function(x) {xmlValue(x)}))
        if (length(abs) == 0) abs = "No Abstract Found"
        
        authorsln = lapply(getNodeSet(efetch, "//LastName"), function(x) {xmlValue(x)});checkAA = unlist(authorsln); authorsini = lapply(getNodeSet(efetch, "//Initials"), function(x) {xmlValue(x)});  if (length(checkAA) == 0) {authorsln = list("No authors name"); authorsini =  list("No   
        authors name")}; authors = unlist(lapply(1:length(authorsln), function(x){return(paste(authorsln[[x]],authorsini[[x]], sep = " "))}))
        arttit = unlist(lapply(getNodeSet(efetch, "//ArticleTitle"), function(x) {xmlValue(x)}))
        if (all == FALSE) return(c(abs, b)) else if (all == TRUE) return(c(arttit,paste(authors, collapse = " "),abs,b))
    })
}
