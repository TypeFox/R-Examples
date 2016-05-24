readFactivaHTML <- FunctionGenerator(function(elem, language, id) {
    function(elem, language, id) {
        encoding <- if(Encoding(elem$content) == "unknown") character(0)
                    else Encoding(elem$content)
        tree <- xmlParse(elem$content, asText=TRUE, encoding=encoding)

        if(is.na(language)) {
            cl <- xmlAttrs(xmlChildren(tree)[[1]])["class"]
            language <- regmatches(cl, regexec("^article ([[:alpha:]]{2})Article$", cl))[[1]][2]
        }

        table <- readHTMLTable(xmlChildren(tree)[[1]])

        # Remove line breaks as paragraphs are used for this
        # (else, line breaks in the source are propagated to the contents)
        text <- gsub("[\n\r]", "",
                     sapply(getNodeSet(tree, "//p[starts-with(@class, 'articleParagraph')]"), xmlValue))
        free(tree)

        # Without this, sometimes table ends up being a mere list
        table <- as.data.frame(table)

        vars <- c("AN", "BY", "CO", "CY", "ED", "HD", "IN", "IPC", "IPD",
                  "LA", "LP", "NS", "PD", "PUB", "RE", "SE", "SN", "TD", "WC")

        # Remove trailing spaces when matching
        data <- as.character(table[match(vars, gsub("[^[A-Z]", "", table[,1])), 2])
        names(data) <- vars

        # Encoding is passed explicitly to work around a bug in XML: htmlParse() and xmlParse() do not set it
        # as they should when asText=TRUE for now
        if(any(Encoding(data) == "unknown"))
            Encoding(data) <- Encoding(elem$content)

        date <- strptime(data[["PD"]], "%d %B %Y")
        if(is.na(date) && isTRUE(data[["PD"]] != "")) {
            # Try C locale, just in case
            old.locale <- Sys.getlocale("LC_TIME")
            Sys.setlocale("LC_TIME", "C")
            date <- strptime(data[["PD"]], "%d %B %Y")
            Sys.setlocale("LC_TIME", old.locale)

            # A bug in Mac OS gives NA when start of month name matches an abbreviated name:
            # http://www.freebsd.org/cgi/query-pr.cgi?pr=141939
            # https://stat.ethz.ch/pipermail/r-sig-mac/2012-June/009296.html
            # Add a workaround for French
            if (Sys.info()["sysname"] == "Darwin")
                date <- strptime(sub("[jJ]uillet", "07", data[["PD"]]), "%d %m %Y")

            if(any(is.na(date)))
                warning(sprintf("Could not parse document date \"%s\". You may need to change the system locale to match that of the corpus. See LC_TIME in ?Sys.setlocale.", data[["PD"]]))
        }

        data[["AN"]] <- gsub("Document ", "", data[["AN"]])

        wc <- as.integer(regmatches(data[["WC"]], regexpr("^[[:digit:]]+", data[["WC"]])))[[1]]

        # Extract useful information: origin, date, and code
        m <- regmatches(data[["AN"]], regexec("^([A-Za-z]+)0*[1-9][0-9]([0-9][0-9][0-3][0-9][0-3][0-9])([A-Za-z0-9])",
                                              data[["AN"]]))[[1]]
        # If extraction failed for some reason, make sure we return a unique identifier
        if(length(m) == 4)
            id <- paste(toupper(m[2]), "-", m[3], "-", m[4], sep="")
        else
            id <- paste(sample(LETTERS, 10), collapse="")

        subject <- if(!is.na(data[["NS"]])) strsplit(data[["NS"]], "( \\| )")[[1]]
                   else character(0)
        # Remove leading code and invisible characters, esp. \n, before matching the pattern
        subject <- gsub("[^[:print:]]", "", subject)
        subject <- gsub(".* : ", "", subject)

        coverage <- if(!is.na(data[["RE"]])) strsplit(data[["RE"]], "( \\| )")[[1]]
                    else character(0)
        # Remove leading code and invisible characters, esp. \n, before matching the pattern
        coverage <- gsub("[^[:print:]]", "", coverage)
        coverage <- gsub(".* : ", "", coverage)

        company <- if(!is.na(data[["CO"]])) strsplit(data[["CO"]], "( \\| )")[[1]]
                   else character(0)
        # Remove leading code and invisible characters, esp. \n, before matching the pattern
        company <- gsub("[^[:print:]]", "", company)
        company <- gsub(".* : ", "", company)

        industry <- if(!is.na(data[["IN"]])) strsplit(data[["IN"]], "( \\| )")[[1]]
                    else character(0)
        # Remove leading code and invisible characters, esp. \n, before matching the pattern
        industry <- gsub("[^[:print:]]", "", industry)
        industry <- gsub(".* : ", "", industry)

        infocode <- if(!is.na(data[["IPC"]])) strsplit(data[["IPC"]], "( \\| )")[[1]]
                    else character(0)
        # Remove leading code and invisible characters, esp. \n, before matching the pattern
        infocode <- gsub("[^[:print:]]", "", infocode)
        infocode <- gsub(".* : ", "", infocode)

        infodesc <- if(!is.na(data[["IPD"]])) strsplit(data[["IPD"]], "( +\\| +| +-+ +| +--+|--+ +|\\._)")[[1]]
                    else character(0)
        # Remove leading code and invisible characters, esp. \n, before matching the pattern
        infodesc <- gsub("[^[:print:]]", "", infodesc)
        infodesc <- gsub(".* : ", "", infodesc)

        # XMLSource uses character(0) rather than NA, do the same
        doc <- PlainTextDocument(x = text,
                                 author = if(!is.na(data[["BY"]])) data[["BY"]] else character(0),
                                 datetimestamp = date,
                                 heading = if(!is.na(data[["HD"]])) data[["HD"]] else character(0),
                                 id = id,
                                 origin = if(!is.na(data[["SN"]])) data[["SN"]] else character(0),
                                 language = language)
        meta(doc, "edition") <- if(!is.na(data[["ED"]])) data[["ED"]] else character(0)
        meta(doc, "section") <- if(!is.na(data[["SE"]])) data[["SE"]] else character(0)
        meta(doc, "subject") <- subject
        meta(doc, "coverage") <- coverage
        meta(doc, "company") <- company
        meta(doc, "industry") <- industry
        meta(doc, "infocode") <- infocode
        meta(doc, "infodesc") <- infodesc
        meta(doc, "wordcount") <- wc
        meta(doc, "publisher") <- if(!is.na(data[["PUB"]])) data[["PUB"]] else character(0)
        meta(doc, "rights") <- if(!is.na(data[["CY"]])) data[["CY"]] else character(0)
        doc
    }
})

