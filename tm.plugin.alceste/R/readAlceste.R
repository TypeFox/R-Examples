readAlceste <- FunctionGenerator(function(elem, language, id) {
    function(elem, language, id) {
        id2 <- regmatches(elem$content[1], regexec("^([[:digit:]]+) \\*", elem$content[1]))[[1]][2]
        # Only override default ID if present
        if(!is.na(id2))
            id <- id2

        starred <- sub("^(\\*\\*\\*\\* +|[[:digit:]]+ \\*)", "", elem$content[1])
        varexpr <- gsub(" ", "", strsplit(starred, "*", starred, fixed=TRUE)[[1]], fixed=TRUE)
        vars <- strsplit(varexpr[nchar(varexpr) > 0], "_", fixed=TRUE)

        # Theme lines, ignored
        skip <- which(grepl("^-\\*", elem$content))

        doc <- PlainTextDocument(x = elem$content[-c(1, skip)],
                                 id = id,
                                 language = language)

        for(v in vars) {
            # Boolean variable (without value after _)
            if(is.na(v[2]))
                meta(doc, v[1]) <- TRUE
            else
                meta(doc, v[1]) <- v[2]
        }

        doc
    }
})

