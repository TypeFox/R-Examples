

isi2bibtex <- function(file) {

    journals <- rbind(c("J. Am. Stat. Assoc.", "Journal of the American Statistical Association", "JASA"),
                      c("J. Stat. Plan. Infer.", "Journal of Statistical Planning and Inference", "JSPI"),
                      c("Biom. J.", "Biometrical Journal", "BJ"),
                      c("Stat. Med.", "Statistics in Medicine", "SiM"))
    colnames(journals) <- c("Abbr", "Title", "ID")

    tfile <- tempfile()
    isitxt <- readLines(file)
    isitxt <- gsub("(^[A-Z][A-Z,0-9])", "\\1:", isitxt, perl = TRUE)
    writeLines(isitxt, con = tfile)

    isidcf <- read.dcf(tfile, fields = c("PT", "AU", "TI", "SO", "LA", "DT", "DE", "ID",
                       "AB", "C1", "RP", "EM", "NR", "TC", "PU", "PI", "PA", "SN",
                       "J9", "JI", "PD", "PY", "VL", "IS", "BP", "EP", "PG", "SC",
                       "GA", "UT"))

    ### journals only
    isidcf <- isidcf[isidcf[,"PT"] == "J",]

    ### missings
    isidcf <- isidcf[!apply(isidcf, 1, function(x) all(is.na(x))),]

    ### rename interesting fields
    cn <- colnames(isidcf)
    colnames(isidcf)[cn == "AU"] <- "author"
    colnames(isidcf)[cn == "TI"] <- "title"
    colnames(isidcf)[cn == "JI"] <- "journal"
    colnames(isidcf)[cn == "PD"] <- "month"
    colnames(isidcf)[cn == "PY"] <- "year"
    colnames(isidcf)[cn == "VL"] <- "volumne"
    colnames(isidcf)[cn == "IS"] <- "number"
    colnames(isidcf)[cn == "UT"] <- "isitag"
    colnames(isidcf)[cn == "DE"] <- "keywords"
    colnames(isidcf)[cn == "TC"] <- "timescited"
    colnames(isidcf)[cn == "AB"] <- "abstract"
    rownames(isidcf) <- 1:nrow(isidcf)
    isidcf[,"title"] <- gsub("\n", " ", isidcf[,"title"])

    ### author names
    for (i in 1:nrow(isidcf)) {
        au <- strsplit(isidcf[i,"author"], "\n")
        names <- strsplit(au[[1]], ", ")
        for (j in 1:length(names))
            names[[j]][2] <- paste(strsplit(names[[j]][2], "")[[1]], ". ", 
                                   sep = "", collapse = "")
        lastnames <- sapply(names, function(x) gsub(" ", "", x[1]))
        if (length(lastnames) > 3) lastnames <- lastnames[1:3]
        jour <- isidcf[i,"journal"]
        indx <- journals[, "Abbr"] == jour
        if (sum(indx) == 1) {
            isidcf[i, "journal"] <- journals[indx, "Title"]
            jkey <- journals[indx, "ID"]
        } else {
            jkey <- gsub("\\.* ", "", jour)
        }
        
        label <- paste(paste(lastnames, collapse = "+"), ":",
                 jkey, ":", isidcf[i,"year"], sep = "")
        rownames(isidcf)[i] <- label
        isidcf[i,"author"] <- paste(sapply(names, function(x) paste(x[2], x[1], sep = "")), 
                                    collapse = " and ")
        title <- isidcf[i, "title"]
        if (!identical(toupper(title), title)) {
            ttmp <- strsplit(title, " ")[[1]]
            lower <- tolower(ttmp) != ttmp
            lower[1] <- FALSE
            ttmp[lower] <- paste("{", ttmp[lower], "}", sep = "")
            isidcf[i, "title"] <- paste(ttmp, collapse = " ")
        }
    }

    tags <- c("author", "title", "journal", "month", "year", "volumne", "number", 
              "isitag", "abstract", "keywords", "timescited")
    isidcf[is.na(isidcf[,"month"]), "month"] <- ""
    for (tag in tags)
        isidcf[,tag] <- paste(tag, " = {", isidcf[,tag], "},", sep = "")
    
    pages <- paste("pages = {", isidcf[, "BP"], "--", isidcf[, "EP"], "},", sep = "")
    headerkey <- paste("@article{", rownames(isidcf), ",", sep = "")

    ret <- vector(mode = "list", length = nrow(isidcf))
    for (i in 1:nrow(isidcf))
        ret[[i]] <- c(headerkey[i], paste("   ", isidcf[i, tags]), 
                      paste("   ", pages[i]), "}", " ")

    unlist(ret)
}
