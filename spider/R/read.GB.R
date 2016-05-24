read.GB <-
function(access.nb, seq.names = access.nb, species.names = TRUE, gene=TRUE, access=TRUE, as.character = FALSE) {
    N <- length(access.nb)
    nrequest <- N%/%400 + as.logical(N%%400)
    X <- character(0)
    for (i in 1:nrequest) {
        a <- (i - 1) * 400 + 1
        b <- 400 * i
        if (i == nrequest) 
            b <- N
        URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=", 
            paste(access.nb[a:b], collapse = ","), "&rettype=gb&retmode=text", 
            sep = "")
        X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
    }
    FI <- grep("^ {0,}ORIGIN", X) + 1
    LA <- which(X == "//") - 1
    obj <- list()
    length(obj) <- N
    for (i in 1:N) {
        tmp <- gsub("[[:digit:] ]", "", X[FI[i]:LA[i]])
        obj[[i]] <- unlist(strsplit(tmp, NULL))
    }
    names(obj) <- seq.names
    if (!as.character) 
        obj <- as.DNAbin(obj)
    if (species.names) {
        tmp <- character(N)
        sp <- grep("ORGANISM", X)
        for (i in 1:N) tmp[i] <- unlist(strsplit(X[sp[i]], " +ORGANISM +"))[2]
        attr(obj, "species") <- gsub(" ", "_", tmp)
    if (gene) {
        tmp2 <- character(N)
        def <- grep("DEFINITION", X)
        for (i in 1:N) tmp2[i] <- unlist(strsplit(X[def[i]], "DEFINITION +"))[2]
        attr(obj, "gene") <- gsub(" ", "_", tmp2)
    if (access) {
        tmp3 <- character(N)
        def <- grep("ACCESSION", X)
        for (i in 1:N) tmp3[i] <- unlist(strsplit(X[def[i]], "ACCESSION +"))[2]
        attr(obj, "accession_num") <- gsub(" ", "_", tmp3)
    }
    names(obj)<-paste(attr(obj,"accession_num"), "|", attr(obj,"species"), sep = "")
    obj
}
}
}

