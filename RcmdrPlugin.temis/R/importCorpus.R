.selectCorpusVariables <- function(source) {
    # Let the user select processing options
    initializeDialog(title=.gettext("Select Variables to Import"))


    vars <- c(.gettext("No variables"), colnames(corpusVars))

    if(source %in% c("factiva", "lexisnexis", "europresse", "twitter"))
        # Keep in sync with import functions
        initialSelection <- which(vars %in% c(.gettext("Origin"), .gettext("Date"), .gettext("Author"),
                                              .gettext("Section"), .gettext("Type"),
                                              .gettext("Time"), .gettext("Truncated"),
                                              .gettext("StatusSource"), .gettext("Retweet"))) - 1
    else
        initialSelection <- seq.int(1, length(vars) - 1)

    varBox <- variableListBox(top, vars,
                              selectmode="multiple",
                              title=.gettext("Select the corpus variables that should be imported:"),
                              initialSelection=initialSelection,
                              listHeight=min(length(vars), 25))

    result <- tclVar()

    onOK <- function() {
        selVars <- getSelection(varBox)

        # In case something goes wrong
        tclvalue(result) <- "error"

        closeDialog()

        # Replace corpusVars with an empty data.frame
        if(length(selVars) == 1 && selVars[1] == .gettext("No variables")) {
            # Because of a bug in Rcmdr, filling the first column with NAs prevents entering data in this columns:
            # use "" instead
            doItAndPrint('corpusVars <- data.frame(var1=factor(rep("", length(corpus))), row.names=names(corpus))')
        }
        # Import only some variables
        else {
            if(length(intersect(vars, selVars)) < length(vars) - 1)
                doItAndPrint(sprintf('corpusVars <- corpusVars[c("%s")]',
                                     paste(setdiff(selVars, .gettext("No variables")),
                                           collapse='", "')))
        }

        tkfocus(CommanderWindow())
        tclvalue(result) <- "success"
    }

    onCancel <- function() {
        if (GrabFocus()) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(CommanderWindow())
        tclvalue(result) <- "cancel"
    }

    OKCancelHelp(helpSubject=importCorpusDlg)
    tkgrid(getFrame(varBox), sticky="nswe", pady=6)
    tkgrid(buttonsFrame, sticky="ew", pady=6)
    # force.wait is required to get checking result to work
    dialogSuffix(force.wait=TRUE)

    return(tclvalue(result) == "success")
}

.langToEncs <- function(lang) {
    switch(lang,
           da="ISO-8859-1, Windows-1252, ISO-8859-15, UTF-8, UTF-16",
           de="ISO-8859-1, Windows-1252, ISO-8859-15, UTF-8, UTF-16",
           en="ASCII, ISO-8859-1, Windows-1252, ISO-8859-15, UTF-8, UTF-16",
           es="ISO-8859-1, Windows-1252, ISO-8859-15, UTF-8, UTF-16",
           fi="ISO-8859-1, Windows-1252, ISO-8859-15, UTF-8, UTF-16",
           fr="ISO-8859-1, Windows-1252, ISO-8859-15, UTF-8, UTF-16",
           hu="ISO-8859-16, ISO-8859-2, Windows-1250, UTF-8, UTF-16",
           it="ISO-8859-1, Windows-1252, ISO-8859-15, UTF-8, UTF-16",
           nl="ISO-8859-1, Windows-1252, ISO-8859-15, UTF-8, UTF-16",
           no="ISO-8859-1, Windows-1252, ISO-8859-15, UTF-8, UTF-16",
           pt="ISO-8859-1, Windows-1252, ISO-8859-15, UTF-8, UTF-16",
           ro="ISO-8859-16, ISO-8859-2, Windows-1250, UTF-8, UTF-16",
           ru="KOI8-R, Windows-1251, ISO-8859-5, UTF-8, UTF-16",
           sv="ISO-8859-1, Windows-1252, ISO-8859-15, UTF-8, UTF-16",
           tr="ISO-8859-9, Windows-1254, UTF-8, UTF-16")
}

# Run all processing steps and extract words list
.processTexts <- function(options, lang) {
        # Check that all texts contain valid characters
        # tolower() is the first function to fail when conversion to lower case is enabled,
        # but if it's not DocumentTermMatrix() calls termFrequencies(), which fails when calling nchar()
        # Since we do not have access to utf8Valid() in R, call nchar() as a way to run a basic validation
        # (a more rigorous check would be to call e.g. grep(), which uses valid_utf8() from PCRE).
        if(inherits(try(for(i in seq(length(corpus))) nchar(corpus[[i]]), silent=TRUE), "try-error")) {
            .Message(sprintf(.gettext("Invalid characters found in document %i. Please check the \"Encoding\" value defined in the import dialog. Most probable encodings for this language are: %s.\n\nIf necessary, use a text editor's \"Save As...\" function to save the corpus in a known encoding."), i, .langToEncs(lang)),
                     type="error")
             return(FALSE)
        }

        doItAndPrint("dtmCorpus <- corpus")

        if(options["twitter"])
            doItAndPrint('dtmCorpus <- tm_map(dtmCorpus, content_transformer(function(x) gsub("http(s?)://[[:alnum:]/\\\\.\\\\-\\\\?=&#_;,]*|\\\\bRT\\\\b", "", x)))')
        if(options["twitter"] && options["removeNames"])
            doItAndPrint('dtmCorpus <- tm_map(dtmCorpus, content_transformer(function(x) gsub("@.+?\\\\b", "", x)))')
        if(options["twitter"] && options["removeHashtags"])
            doItAndPrint('dtmCorpus <- tm_map(dtmCorpus, content_transformer(function(x) gsub("#.+?\\\\b", "", x)))')

        if(options["lowercase"])
            doItAndPrint("dtmCorpus <- tm_map(dtmCorpus, content_transformer(tolower))")

        if(options["punctuation"]) {
            # The default tokenizer does not get rid of punctuation *and of line breaks!*, which
            # get concatenated with surrounding words
            # This also avoids French articles and dash-linked words from getting concatenated with their noun
            doItAndPrint("dtmCorpus <- tm_map(dtmCorpus, content_transformer(function(x) gsub(\"([\'\U2019\\n\U202F\U2009]|[[:punct:]]|[[:space:]]|[[:cntrl:]])+\", \" \", x)))")
        }

        if(options["digits"])
            doItAndPrint("dtmCorpus <- tm_map(dtmCorpus, removeNumbers)")

        return(TRUE)
}

.buildDictionary <- function(stemming, customStemming, lang) {
    if(stemming) {
        doItAndPrint("library(SnowballC)")
        doItAndPrint(sprintf('dictionary <- data.frame(row.names=colnames(dtm), "%s"=col_sums(dtm), "%s"=wordStem(colnames(dtm), "%s"), "%s"=ifelse(colnames(dtm) %%in%% stopwords("%s"), "%s", ""), stringsAsFactors=FALSE)',
                             .gettext("Occurrences"), .gettext("Stemmed.Term"), lang, .gettext("Stopword"), lang, .gettext("Stopword")))
    }
    else if(customStemming){
        doItAndPrint(sprintf('dictionary <- data.frame(row.names=colnames(dtm), "%s"=col_sums(dtm), "%s"=colnames(dtm), "%s"=ifelse(colnames(dtm) %%in%% stopwords("%s"), "%s", ""), stringsAsFactors=FALSE)',
                             .gettext("Occurrences"), .gettext("Stemmed.Term"), .gettext("Stopword"), lang, .gettext("Stopword")))
    }
    else {
        doItAndPrint(sprintf('dictionary <- data.frame(row.names=colnames(dtm), "%s"=col_sums(dtm), "%s"=ifelse(colnames(dtm) %%in%% stopwords("%s"), "%s", ""), stringsAsFactors=FALSE)',
                             .gettext("Occurrences"), .gettext("Stopword"), lang, .gettext("Stopword")))
    }
}

.prepareDtm <- function(stopwords, stemming, customStemming, lang) {
    if(customStemming) {
        if(stopwords) doItAndPrint(sprintf('dictionary[Terms(dtm) %%in%% stopwords("%s"), "%s"] <- ""',
                                           lang, .gettext("Stemmed.Term")))

        doItAndPrint("dictionary <- editDictionary(dictionary)")
        doItAndPrint(sprintf('dtm <- rollup(dtm, 2, dictionary[["%s"]])', .gettext("Stemmed.Term")))
        doItAndPrint('dtm <- dtm[, Terms(dtm) != ""]')
    }
    else if(stemming && stopwords) {
        doItAndPrint(sprintf('dtm <- dtm[, !colnames(dtm) %%in%% stopwords("%s")]', lang))
        doItAndPrint("dtm <- rollup(dtm, 2, dictionary[colnames(dtm), 2])")
    }
    else if(stemming) {
        doItAndPrint("dtm <- rollup(dtm, 2, dictionary[[2]])")
    }
    else if(stopwords) {
        doItAndPrint(sprintf('dtm <- dtm[, !colnames(dtm) %%in%% stopwords("%s")]', lang))
    }

    doItAndPrint('attr(dtm, "dictionary") <- dictionary')
    doItAndPrint("rm(dictionary)")
}

importCorpusDlg <- function() {
    # Let the user select processing options
    initializeDialog(title=.gettext("Import Corpus"))

    setState <- function(...) {
        if(tclvalue(sourceVariable) %in% c("dir", "file", "alceste")) {
            tkconfigure(comboEnc, state="normal")

            if(tclvalue(tclEnc) == "UTF-8")
                tclvalue(tclEnc) <- autoEnc
        }
        else {
            tkconfigure(comboEnc, state="disabled")

            if(tclvalue(tclEnc) == autoEnc)
                tclvalue(tclEnc) <- "UTF-8"
        }
    }

    radioButtons(name="source",
                 buttons=c("dir", "file", "factiva", "lexisnexis", "europresse", "alceste", "twitter"),
                 labels=c(.gettext("Directory containing plain text files"),
                          .gettext("Spreadsheet file (CSV, XLS, ODS...)"),
                          .gettext("Factiva XML or HTML file(s)"),
                          .gettext("LexisNexis HTML file(s)"),
                          .gettext("Europresse HTML file(s)"),
                          .gettext("Alceste file(s)"),
                          .gettext("Twitter search")),
                 title=.gettext("Load corpus from:"),
                 right.buttons=FALSE,
                 command=setState)

    autoEnc <- .gettext("detect automatically")
    tclEnc <- tclVar(autoEnc)
    # Do not use state="readonly" since it may be easier to type the encoding name by hand
    # than choose it in the long list
    comboEnc <- ttkcombobox(top, width=20, textvariable=tclEnc,
                            values=c(autoEnc, iconvlist()))

    # Keep in sync with .processTexts()
    # TRANSLATORS: replace 'en' with your language's ISO 639 two-letter code
    languages <- c(da="Dansk (da)", de="Deutsch (de)", en="English (en)", es="Espa\u00F1ol (es)",
                   fi="Suomi (fi)", fr="Fran\u00E7ais (fr)", hu="Magyar (hu)", it="Italiano (it)",
                   nl="Nederlands (nl)", no="Norsk (no)", pt="Portugu\u0EAs (pt)",
                   ro="Rom\u00E2n\u0103 (ro)",
                   ru="\u0440\u0443\u0441\u0441\u043A\u0438\u0439 \u044F\u0437\u044B\u043A (ru)",
                   sv="Svenska (sv)", tr="T\u00FCrk\u00E7e (tr)")
    tclLang <- tclVar(languages[.gettext("en")])
    comboLang <- ttkcombobox(top, width=20, textvariable=tclLang, state="readonly", values=languages)

    tk2tip(comboEnc, sprintf(.gettext("Most probable encodings for this language:\n%s"),
                             .langToEncs(.gettext("en"))))
    tkbind(comboLang, "<<ComboboxSelected>>", function() {
        # On Windows when the current locale does not support characters in the language name,
        # the value retrieved using tclvalue() contains invalid characters
        # The only way to find the language is to extract the relevant two ASCII letters
        rawLang <- tclvalue(tclLang)
        lang <- substring(rawLang, nchar(rawLang) - 2, nchar(rawLang) - 1)
        tk2tip(comboEnc, sprintf(.gettext("Most probable encodings for this language:\n%s"),
                                 .langToEncs(lang)))
    })

    checkBoxes(frame="processingFrame",
               boxes=c("lowercase", "punctuation", "digits", "stopwords",
                       "stemming", "customStemming"),
               initialValues=c(1, 1, 1, 0, 1, 0),
               # Keep in sync with strings in initOutputFile()
               labels=c(.gettext("Ignore case"), .gettext("Remove punctuation"),
                        .gettext("Remove digits"), .gettext("Remove stopwords"),
                        .gettext("Apply stemming"), .gettext("Edit stemming manually")),
               title=.gettext("Text processing:"))

    tclChunks <- tclVar(0)
    tclNParagraphs <- tclVar(1)
    chunksButton <- ttkcheckbutton(top, variable=tclChunks,
                                   text=.gettext("Split texts into smaller documents"),
                                    command=function() {
                                        if(tclvalue(tclChunks) == 1)
                                            tkconfigure(chunksSlider, state="active")
                                        else
                                            tkconfigure(chunksSlider, state="disabled")
                                    })
    chunksSlider <- tkscale(top, from=1, to=20, showvalue=TRUE, variable=tclNParagraphs,
                            resolution=1, orient="horizontal", state="disabled")


    onOK <- function() {
        source <- tclvalue(sourceVariable)
        twitter <- source == "twitter"
        lowercase <- tclvalue(lowercaseVariable) == 1
        punctuation <- tclvalue(punctuationVariable) == 1
        digits <- tclvalue(digitsVariable) == 1
        stopwords <- tclvalue(stopwordsVariable) == 1
        stemming <- tclvalue(stemmingVariable) == 1
        customStemming <- tclvalue(customStemmingVariable) == 1

        # On Windows when the current locale does not support characters in the language name,
        # the value retrieved using tclvalue() contains invalid characters
        # The only way to find the language is to extract the relevant two ASCII letters
        rawLang <- tclvalue(tclLang)
        lang <- substring(rawLang, nchar(rawLang) - 2, nchar(rawLang) - 1)

        enc <- tclvalue(tclEnc)
        if(enc == autoEnc) enc <- ""

        if(enc != "" && !enc %in% iconvlist()) {
            .Message(.gettext('Unsupported encoding: please select an encoding from the list.'),
                     "error", parent=top)
            return()
        }

        closeDialog()

        setBusyCursor()
        on.exit(setIdleCursor())

	if(stemming) {
            # If we do not close the dialog first, the CRAN mirror chooser will not respond
	    if(!.checkAndInstall("SnowballC", .gettext("Package SnowballC is needed to perform stemming. Do you want to install it?\n\nThis requires a working Internet connection.")))
                return()
        }

        # Remove objects left from a previous analysis to avoid confusion
        # (we assume later existing objects match the current corpus)
        objects <- c("corpus", "corpusVars", "dtm", "wordsDtm", "lengthsDtm", "voc", "lengths", "coocs",
                     "termFreqs", "absTermFreqs", "varTermFreqs", "freqTerms", "specTerms", "docSeries",
                     ".last.table", "corpusClust", "corpusSubClust", "corpusCa", "plottingCa")
        if(any(sapply(objects, exists))) {
            # Needed to avoid a warning about corpusVars not being available
            # This should be removed once we can depend on the new version of Rcmdr accepting ActiveDataSet(NULL)
            if(exists("corpusVars")) {
                putRcmdr(".activeDataSet", NULL)
		        Variables(NULL)
		        Numeric(NULL)
		        Factors(NULL)
		        RcmdrTclSet("dataSetName", gettextRcmdr("<No active dataset>"))
		        putRcmdr(".activeModel", NULL)
		        RcmdrTclSet("modelName", gettextRcmdr("<No active model>"))
		        tkconfigure(getRcmdr("dataSetLabel"), foreground="red")
		        tkconfigure(getRcmdr("modelLabel"), foreground="red")
            }

            doItAndPrint(paste("rm(", paste(objects[sapply(objects, exists)], collapse=", "), ")", sep=""))
            gc()
        }

        HTMLSetFile(NULL)
        activateMenus()

        # Import corpus
        res <- switch(source,
                      dir=importCorpusFromDir(lang, enc),
                      file=importCorpusFromFile(lang, enc),
                      factiva=importCorpusFromFactiva(lang),
                      lexisnexis=importCorpusFromLexisNexis(lang),
                      europresse=importCorpusFromEuropresse(lang),
                      alceste=importCorpusFromAlceste(lang, enc),
                      twitter=importCorpusFromTwitter(lang))

        # If loading failed, do not add errors to errors
        if(!(isTRUE(res) || is.list(res)) || length(corpus) == 0)
            return()

        # If source-specific functions load variables, they create corpusVars; else, create an empty data frame
        if(exists("corpusVars")) {
            if(!.selectCorpusVariables(source)) return()
        }
        else {
            # Because of a bug in Rcmdr, filling the first column with NAs prevents entering data in this columns:
            # use "" instead
            doItAndPrint('corpusVars <- data.frame(var1=factor(rep("", length(corpus))), row.names=names(corpus))')
        }

        # Needed because functions above set it to idle on exit
        setBusyCursor()
        on.exit(setIdleCursor())

        doItAndPrint('activeDataSet("corpusVars")')
        doItAndPrint("setCorpusVariables()")

        # Create chunks
        if(tclvalue(tclChunks) == 1) {
            doItAndPrint(sprintf("corpus <- splitTexts(corpus, %s)", tclvalue(tclNParagraphs)))
            doItAndPrint('meta(corpus, type="corpus", tag="split") <- TRUE')
        }

        # Process texts
        if(!.processTexts(c(twitter=twitter, lowercase=lowercase, punctuation=punctuation,
                            digits=digits, stopwords=stopwords,
                            stemming=stemming, customStemming=customStemming,
                            removeHashtags=res$removeHashtags, removeNames=res$removeNames),
                          lang))
            return()

        doItAndPrint("dtm <- DocumentTermMatrix(dtmCorpus, control=list(tolower=FALSE, wordLengths=c(2, Inf)))")
        doItAndPrint("rm(dtmCorpus)")

        if(!exists("dtm") || is.null(dtm)) {
            return()
        }
        else if(nrow(dtm) == 0 || ncol(dtm) == 0) {
            doItAndPrint("dtm")
            return()
        }

        .buildDictionary(stemming, customStemming, lang)
        .prepareDtm(stopwords, stemming, customStemming, lang)

        gc()


        # Language is used again when creating the dtm to analyse word lengths
        doItAndPrint(sprintf('meta(corpus, type="corpus", tag="language") <- attr(dtm, "language") <- "%s"', lang))

        # Do not call doItAndPrint() until the bug in Rcmdr is fixed (long quoted strings create problems
        # when splitting commands when more than one pair of quotes is present)
        justDoIt(sprintf('meta(corpus, type="corpus", tag="source") <- "%s"', res$source))

        doItAndPrint(sprintf('meta(corpus, type="corpus", tag="processing") <- attr(dtm, "processing") <- c(lowercase=%s, punctuation=%s, digits=%s, stopwords=%s, stemming=%s, customStemming=%s, twitter=%s, removeHashtags=%s, removeNames=%s)',
                             lowercase, punctuation,
                             digits, stopwords, stemming, customStemming, twitter,
                             ifelse(is.null(res$removeHashtags), NA, res$removeHashtags),
                             ifelse(is.null(res$removeNames), NA, res$removeNames)))

        doItAndPrint("corpus")
        doItAndPrint("dtm")

        activateMenus()

        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="importCorpusDlg")
    tkgrid(sourceFrame, columnspan=3, sticky="w", pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Language of texts in the corpus:")), sticky="w", pady=6)
    tkgrid(comboLang, sticky="ew", pady=6, row=1, column=1, columnspan=2)
    tkgrid(labelRcmdr(top, text=.gettext("File encoding:")), sticky="w", pady=6)
    tkgrid(comboEnc, sticky="ew", pady=6, row=2, column=1, columnspan=2)
    tkgrid(.titleLabel(top, text=.gettext("Text splitting:")),
           sticky="ws", pady=c(12, 6))
    tkgrid(chunksButton, sticky="w", pady=6, columnspan=3)
    tkgrid(labelRcmdr(top, text=.gettext("Size of new documents:")),
           chunksSlider, labelRcmdr(top, text=.gettext("paragraphs")), sticky="w", pady=6)
    tkgrid(processingFrame, columnspan=3, sticky="w", pady=6)
    tkgrid(buttonsFrame, columnspan=3, sticky="ew", pady=6)
    tkgrid.columnconfigure(top, 0, pad=12)
    tkgrid.columnconfigure(top, 1, pad=12)
    tkgrid.columnconfigure(top, 2, pad=12)
    dialogSuffix(focus=comboLang)
}

# Choose a directory to load texts from
importCorpusFromDir <- function(language=NA, encoding="") {
    dir <- tclvalue(tkchooseDirectory(initialdir=getwd(),
                                      parent=CommanderWindow()))
    if (dir == "") return(FALSE)

    setBusyCursor()
    on.exit(setIdleCursor())

    if(!is.na(language))
        language <- paste("\"", language, "\"", sep="")

    if(encoding == "") {
        encs <- table(sapply(list.files(dir, full.names=TRUE),
                             function(f) stri_enc_detect(readBin(f, "raw", 50000))[[1]]$Encoding[1]))
        encoding <- names(encs)[order(encs, decreasing=TRUE)][1]
    }

    if(is.null(encoding))
        encoding <- ""

    doItAndPrint(sprintf('corpus <- Corpus(DirSource("%s", encoding="%s"), readerControl=list(language=%s))',
                         dir, encoding, language))

    list(source=sprintf(.gettext("directory %s"), dir))
}

# Choose a CSV file to load texts and variables from
importCorpusFromFile <- function(language=NA, encoding="") {
    file <- tclvalue(tkgetOpenFile(filetypes=sprintf("{{%s} {.csv .CSV .tsv .TSV .dbf .DBF .ods .ODS .xls .XLS .xlsx .XLSX .mdb .MDB .accdb .ACCDB}} {{%s} {.csv .CSV}} {{%s} {.tsv .txt .dat .TSV .TXT .DAT}} {{%s} {.dbf .DBF}} {{%s} {.ods .ODS}} {{%s} {.xls .XLS}} {{%s} {.xlsx .XLSX}} {{%s} {.mdb .MDB}} {{%s} {.accdb .ACCDB}} {{%s} {*}}",
                                                     .gettext("All supported types"),
                                                     .gettext("Comma-separated values (CSV) file"),
                                                     .gettext("Tab-separated values (TSV) file"),
                                                     .gettext("dBase file"),
                                                     .gettext("Open Document Spreadsheet file"),
                                                     .gettext("Excel file"),
                                                     .gettext("Excel 2007 file"),
                                                     .gettext("Access database"),
                                                     .gettext("Access 2007 database"),
                                                     .gettext("All files")),
                                   parent=CommanderWindow()))

    if (file == "") return(FALSE)

    setBusyCursor()
    on.exit(setIdleCursor())

    # Code adapted from Rcommander's data-menu.R file
    # The following function was contributed by Matthieu Lesnoff
    # (added with small changes by J. Fox, 20 July 06 & 30 July 08)
    # Licensed under GNU GPL (version >= 2)
    sop <- match(".", rev(strsplit(file, NULL)[[1]]))[1]
    ext <- tolower(substring(file, nchar(file) - sop + 2, nchar(file)))

    if(ext == "csv") {
        # Try to guess the separator from the most common character of ; and ,
        # This should work in all cases where text is not too long
        excerpt <- readLines(file, 50)
        n1 <- sum(sapply(gregexpr(",", excerpt), length))
        n2 <- sum(sapply(gregexpr(";", excerpt), length))

        if(encoding == "")
            encoding <- stri_enc_detect(readBin(file, "raw", 1024))[[1]]$Encoding[1]

        if(is.null(encoding))
            encoding <- ""

        if(n1 > n2)
            doItAndPrint(sprintf('corpusDataset <- read.csv("%s", fileEncoding="%s")', file, encoding))
        else
            doItAndPrint(sprintf('corpusDataset <- read.csv2("%s", fileEncoding="%s")', file, encoding))
    }
    else if(ext %in% c("tsv", "txt", "dat")) {
        doItAndPrint(sprintf('corpusDataset <- read.delim("%s", fileEncoding="%s")', file, encoding))
    }
    else if(ext == "dbf") {
        Library(foreign)
        doItAndPrint(paste("corpusDataset <- read.dbf(\"", file, "\")", sep=""))
    }
    else if(ext == "ods") {
        # ROpenOffice is not available as binary, thus most likely to fail on Windows and Mac OS
        if(!"ROpenOffice" %in% rownames(available.packages(contrib.url("http://www.omegahat.org/R/")))) {
	    .Message(.gettext("Loading OpenDocument spreadsheets (.ods) is not supported on your system.\nYou should save your data set as a CSV file or as an Excel spreadsheet (.xls)."),
                 type="error")
            return(FALSE)
        }
	    else if(!requireNamespace("ROpenOffice")) {
            response <- tkmessageBox(message=.gettext("Loading OpenDocument spreadsheets (.ods) requires the ROpenOffice package.\nDo you want to install it?"),
                                     icon="question", type="yesno")

            if (tclvalue(response) == "yes")
	            install.packages("ROpenOffice", repos="http://www.omegahat.org/R", type="source")
            else
                return(FALSE)
        }

        doItAndPrint("library(ROpenOffice)")
        doItAndPrint(paste("corpusDataset <- read.ods(\"", file, "\")", sep=""))
    }
    else if(ext %in% c("xls", "xlsx", "mdb", "accdb")) {
        if(.Platform$OS.type != "windows") {
	    .Message(.gettext("Loading Excel and Access files is only supported on Windows.\nYou should save your data set as a CSV file or as an OpenDocument spreadsheet (.ods)."),
                 type="error")
            return(FALSE)
        }
	else if(!.checkAndInstall("RODBC", .gettext("The RODBC package is needed to read Excel and Access files.\nDo you want to install it?"))) {
            return(FALSE)
        }
        else if(!any(grepl(ext, odbcDataSources()))) {
	    .Message(.gettext("No ODBC driver for this file type was found.\nYou probably need to install Excel or Access, or separate ODBC drivers."),
                 type="error")
            return(FALSE)
        }

        doItAndPrint("library(RODBC)")

        channelStr <- switch(EXPR = ext,
        	             xls = "odbcConnectExcel",
        	             xlsx = "odbcConnectExcel2007",
        	             mdb = "odbcConnectAccess",
        	             accdb = "odbcConnectAccess2007")
        doItAndPrint(paste("channel <- ", channelStr, "(\"", file, "\")", sep=""))

        # For Excel and Access, need to select a particular sheet or table
        tabdat <- sqlTables(channel)
        names(tabdat) <- tolower(names(tabdat))

        if(ext == "mdb" || ext == "accdb")
            tabdat <- tabdat[tabdat$table_type == "TABLE", 3]

        if(ext == "xls" || ext == "xlsx") {
            tabname <- tabdat$table_name
            tabdat <- ifelse(tabdat$table_type == "TABLE",
                             substring(tabname, 2, nchar(tabname) - 2),
                             substring(tabname, 1, nchar(tabname) - 1))
        }

        # If there are several tables
        if(length(tabdat) > 1)
            fil <- tk_select.list(sort(tabdat), title=.gettext("Select one table"))
        else
            fil <- tabdat

        if(fil == "") {
            .Message(.gettext("No table selected"), type="error")
            return(FALSE)
        }

        if(ext == "xls" || ext == "xlsx")
            fil <- paste("[", fil, "$]", sep = "")

        # Retrieve the data
        command <- sprintf('sqlQuery(channel=channel, "select * from %s")', fil)
        doItAndPrint(paste("corpusDataset <- ", command, sep = ""))
        doItAndPrint("odbcCloseAll()")
    }
    else {
        .Message(.gettext("File has unknown extension, assuming it is in the tab-separated values format."), "warning")
        doItAndPrint(paste("corpusDataset <- read.delim(\"", file, "\")", sep=""))
    }

    # In case something went wrong, no point in continuing
    if(is.null(corpusDataset))
        return(FALSE)

    initializeDialog(title=.gettext("Select Text Variable"))
    varBox <- variableListBox(top, names(corpusDataset),
                              selectmode="single",
                              initialSelection=0,
                              title=.gettext("Select the variable containing the text of documents"))

    var <- NULL

    onOK <- function() {
        var <<- getSelection(varBox)

        closeDialog()
        tkfocus(CommanderWindow())
    }

    onCancel <- function() {
        closeDialog()
        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject=importCorpusDlg)
    tkgrid(getFrame(varBox), sticky="nswe", pady=6)
    tkgrid(buttonsFrame, sticky="ew", pady=6)
    dialogSuffix(force.wait=TRUE)

    if(is.null(var))
        return(FALSE)

    if(!is.na(language))
        language <- paste("\"", language, "\"", sep="")

    doItAndPrint(sprintf('corpus <- Corpus(DataframeSource(corpusDataset["%s"]), readerControl=list(language=%s))',
                         var, language))

    if(ncol(corpusDataset) > 1)
        doItAndPrint(sprintf('corpusVars <- corpusDataset[!names(corpusDataset) == "%s"]', var))

    list(source=sprintf(.gettext("spreadsheet file %s"), file))
}

# Extract local per-document meta-data and return a data frame
extractMetadata <- function(corpus, date=TRUE) {
    if(date) {
        dates <- lapply(corpus, meta, "datetimestamp")
        dates <- sapply(dates, function(x) if(length(x) > 0) as.character(x) else NA)
        vars <- data.frame(origin=rep(NA, length(corpus)),
                           date=dates,
                           author=rep(NA, length(corpus)),
                           section=rep(NA, length(corpus)))
    }
    else {
        vars <- data.frame(origin=rep(NA, length(corpus)),
                           author=rep(NA, length(corpus)),
                           section=rep(NA, length(corpus)))
    }

    specialTags <- c("subject", "coverage", "company", "stocksymbol", "industry", "infocode", "infodesc")

    tags <- setdiff(unique(unlist(lapply(corpus, function(x) names(meta(x))))),
                    c("datetimestamp", "heading", "id", "language", specialTags))
    for(tag in tags) {
        var <- lapply(corpus, meta, tag)
        # paste() is here to prevent an error in case x contains more than one elemen
        # This typically happens with Rights
        var <- lapply(var, function(x) if(length(x) > 0) paste(x, collapse=" ") else NA)
        vars[[tag]] <- unlist(var)
    }

    # Keep in sync with importCorpusFromTwitter()
    colnames(vars)[colnames(vars) == "origin"] <- .gettext("Origin")
    colnames(vars)[colnames(vars) == "date"] <- .gettext("Date")
    colnames(vars)[colnames(vars) == "author"] <- .gettext("Author")
    colnames(vars)[colnames(vars) == "section"] <- .gettext("Section")
    colnames(vars)[colnames(vars) == "type"] <- .gettext("Type")
    colnames(vars)[colnames(vars) == "edition"] <- .gettext("Edition")
    colnames(vars)[colnames(vars) == "wordcount"] <- .gettext("Word.Count")
    colnames(vars)[colnames(vars) == "pages"] <- .gettext("Pages")
    colnames(vars)[colnames(vars) == "publisher"] <- .gettext("Publisher")
    colnames(vars)[colnames(vars) == "rights"] <- .gettext("Rights")

    # Drop variables with only NAs, which can appear with sources that do not support them
    vars <- vars[sapply(vars, function(x) sum(!is.na(x))) > 0]


    # Tags that contain several values and have to be represented using dummies
    meta <- sapply(corpus, function(x) meta(x)[specialTags])
    # Tags missing from all documents
    meta <- meta[!is.na(rownames(meta)),]
    # Tags missing from some documents
    meta[] <- sapply(meta, function(x) if(is.null(x)) NA else x)

    for(tag in rownames(meta)) {
        var <- meta[tag,]
        levs <- unique(unlist(var))
        levs <- levs[!is.na(levs)]

        if(length(levs) == 0)
            next

        # We remove the identifier before ":" and abbreviate the names since they can get out of control
        for(lev in levs)
            vars[[make.names(substr(gsub("^[[:alnum:]]+ : ", "", lev), 1, 20))]] <- sapply(var, function(x) lev %in% x)
    }

    rownames(vars) <- names(corpus)

    vars
}

# Choose a Factiva XML or HTML file to load texts and variables from
importCorpusFromFactiva <- function(language=NA) {
    if(!.checkAndInstall("tm.plugin.factiva",
                         .gettext("The tm.plugin.factiva package is needed to import corpora from Factiva files.\nDo you want to install it?")))
        return(FALSE)

    filestr <- tclvalue(tkgetOpenFile(filetypes=sprintf("{{%s} {.xml .htm .html .aspx .XML .HTM .HTML .ASPX}} {{%s} {*}}",
                                                        .gettext("Factiva XML and HTML files"),
                                                        .gettext("All files")),
                                      multiple=TRUE,
                                      parent=CommanderWindow()))

    if (filestr == "") return(FALSE)

    setBusyCursor()
    on.exit(setIdleCursor())

    doItAndPrint("library(tm.plugin.factiva)")

    # tkgetOpenFile() is terrible: if path contains a space, file paths are surrounded by {}
    # If no spaces are present, they are not, but in both cases the separator is a space
    if(substr(filestr, 0, 1) == "{")
        files <- gsub("\\{|\\}", "", strsplit(filestr, "\\} \\{")[[1]])
    else
        files <- strsplit(filestr, " ")[[1]]

    if(!is.na(language))
        language <- paste("\"", language, "\"", sep="")

    if(length(files) == 1) {
        doItAndPrint(sprintf("corpus <- Corpus(FactivaSource(\"%s\"), readerControl=list(language=%s))",
                             files[1], language))

		if(!exists("corpus") || length(corpus) == 0) {
		    .Message(.gettext("Reading the specified file failed. Are you sure this file is in the correct format?"),
		             type="error")

		    return(FALSE)
		}
    }
    else {
        doItAndPrint(sprintf('corpusList <- vector(%s, mode="list")', length(files)))
        for(i in seq(length(files))) {
            doItAndPrint(sprintf('corpusList[[%s]] <- Corpus(FactivaSource("%s"), readerControl=list(language=%s))',
                                 i, files[i], language))

			if(length(corpusList[[i]]) == 0) {
				.Message(.gettext("Reading the specified file failed. Are you sure this file is in the correct format?"),
				         type="error")

				return(FALSE)
			}
        }

        doItAndPrint("corpus <- do.call(c, c(corpusList, list(recursive=TRUE)))")
        doItAndPrint("rm(corpusList)")
        gc()
    }

    # We rely on names/IDs later e.g. in showCorpusCa() because we cannot use indexes when documents are skipped
    # In rare cases, duplicated IDs can happen since Factiva plugin truncates them: ensure they are unique
    doItAndPrint("names(corpus) <- make.unique(names(corpus))")

    doItAndPrint("corpusVars <- extractMetadata(corpus)")

    list(source=sprintf(.ngettext(length(files), "Factiva file %s", "Factiva files %s"),
                        paste(files, collapse=", ")))
}

# Choose a LexisNexis HTML file to load texts and variables from
importCorpusFromLexisNexis <- function(language=NA) {
    if(!.checkAndInstall("tm.plugin.lexisnexis",
                         .gettext("The tm.plugin.lexisnexis package is needed to import corpora from Factiva files.\nDo you want to install it?")))
        return(FALSE)

    filestr <- tclvalue(tkgetOpenFile(filetypes=sprintf("{{%s} {.xml .htm .html .aspx .XML .HTM .HTML .ASPX}} {{%s} {*}}",
                                                        .gettext("LexisNexis HTML files"),
                                                        .gettext("All files")),
                                      multiple=TRUE,
                                      parent=CommanderWindow()))

    if (filestr == "") return(FALSE)

    setBusyCursor()
    on.exit(setIdleCursor())

    doItAndPrint("library(tm.plugin.lexisnexis)")

    # tkgetOpenFile() is terrible: if path contains a space, file paths are surrounded by {}
    # If no spaces are present, they are not, but in both cases the separator is a space
    if(substr(filestr, 0, 1) == "{")
        files <- gsub("\\{|\\}", "", strsplit(filestr, "\\} \\{")[[1]])
    else
        files <- strsplit(filestr, " ")[[1]]

    if(!is.na(language))
        language <- paste("\"", language, "\"", sep="")

    if(length(files) == 1) {
        doItAndPrint(sprintf("corpus <- Corpus(LexisNexisSource(\"%s\"), readerControl=list(language=%s))",
                             files[1], language))

		if(!exists("corpus") || length(corpus) == 0) {
		    .Message(.gettext("Reading the specified file failed. Are you sure this file is in the correct format?"),
		             type="error")

		    return(FALSE)
		}
    }
    else {
        doItAndPrint(sprintf('corpusList <- vector(%s, mode="list")', length(files)))
        for(i in seq(length(files))) {
            doItAndPrint(sprintf('corpusList[[%s]] <- Corpus(LexisNexisSource("%s"), readerControl=list(language=%s))',
                                 i, files[i], language))

			if(length(corpusList[[i]]) == 0) {
				.Message(.gettext("Reading the specified file failed. Are you sure this file is in the correct format?"),
				         type="error")

				return(FALSE)
			}
        }

        doItAndPrint("corpus <- do.call(c, c(corpusList, list(recursive=TRUE)))")
        doItAndPrint("rm(corpusList)")
        gc()
    }

    # We rely names/IDs this later e.g. in showCorpusCa() because we cannot use indexes when documents are skipped
    # In rare cases, duplicated IDs can happen since LexisNexis does not provide any identifier
    # this is unlikely, though, since we include in the ID the document number in the corpus
    doItAndPrint("names(corpus) <- make.unique(names(corpus))")

    doItAndPrint("corpusVars <- extractMetadata(corpus)")

    list(source=sprintf(.ngettext(length(files), "LexisNexis file %s", "LexisNexis files %s"),
                        paste(files, collapse=", ")))
}

# Choose a Europresse HTML file to load texts and variables from
importCorpusFromEuropresse <- function(language=NA) {
    if(!.checkAndInstall("tm.plugin.europresse",
                         .gettext("The tm.plugin.europresse package is needed to import corpora from Europresse files.\nDo you want to install it?")))
        return(FALSE)

    filestr <- tclvalue(tkgetOpenFile(filetypes=sprintf("{{%s} {.htm .html .aspx .HTM .HTML .ASPX}} {{%s} {*}}",
                                                        .gettext("Europresse HTML files"),
                                                        .gettext("All files")),
                                      multiple=TRUE,
                                      parent=CommanderWindow()))

    if (filestr == "") return(FALSE)

    setBusyCursor()
    on.exit(setIdleCursor())

    doItAndPrint("library(tm.plugin.europresse)")

    # tkgetOpenFile() is terrible: if path contains a space, file paths are surrounded by {}
    # If no spaces are present, they are not, but in both cases the separator is a space
    if(substr(filestr, 0, 1) == "{")
        files <- gsub("\\{|\\}", "", strsplit(filestr, "\\} \\{")[[1]])
    else
        files <- strsplit(filestr, " ")[[1]]

    if(!is.na(language))
        language <- paste("\"", language, "\"", sep="")

    if(length(files) == 1) {
        doItAndPrint(sprintf("corpus <- Corpus(EuropresseSource(\"%s\"), readerControl=list(language=%s))",
                             files[1], language))

		if(!exists("corpus") || length(corpus) == 0) {
		    .Message(.gettext("Reading the specified file failed. Are you sure this file is in the correct format?"),
		             type="error")

		    return(FALSE)
		}
    }
    else {
        doItAndPrint(sprintf('corpusList <- vector(%s, mode="list")', length(files)))
        for(i in seq(length(files))) {
            doItAndPrint(sprintf('corpusList[[%s]] <- Corpus(EuropresseSource("%s"), readerControl=list(language=%s))',
                                 i, files[i], language))

			if(length(corpusList[[i]]) == 0) {
				.Message(.gettext("Reading the specified file failed. Are you sure this file is in the correct format?"),
				         type="error")

				return(FALSE)
			}
        }

        doItAndPrint("corpus <- do.call(c, c(corpusList, list(recursive=TRUE)))")
        doItAndPrint("rm(corpusList)")
        gc()
    }

    # We rely names/IDs this later e.g. in showCorpusCa() because we cannot use indexes when documents are skipped
    # In rare cases, duplicated IDs can happen since Europresse plugin truncates them: ensure they are unique
    doItAndPrint("names(corpus) <- make.unique(substr(names(corpus), 1, 20))")

    doItAndPrint("corpusVars <- extractMetadata(corpus)")

    list(source=sprintf(.ngettext(length(files), "Europresse file %s", "Europresse files %s"),
                        paste(files, collapse=", ")))
}


# Choose an Alceste file to load texts and variables from
importCorpusFromAlceste <- function(language=NA, encoding="") {
    if(!.checkAndInstall("tm.plugin.alceste",
                         .gettext("The tm.plugin.alceste package is needed to import corpora from Alceste files.\nDo you want to install it?")))
        return(FALSE)

    filestr <- tclvalue(tkgetOpenFile(filetypes=sprintf("{{%s} {.txt .TXT}} {{%s} {*}}",
                                                        .gettext("Alceste files"),
                                                        .gettext("All files")),
                                      multiple=TRUE,
                                      parent=CommanderWindow()))

    if (filestr == "") return(FALSE)

    setBusyCursor()
    on.exit(setIdleCursor())

    doItAndPrint("library(tm.plugin.alceste)")

    # tkgetOpenFile() is terrible: if path contains a space, file paths are surrounded by {}
    # If no spaces are present, they are not, but in both cases the separator is a space
    if(substr(filestr, 0, 1) == "{")
        files <- gsub("\\{|\\}", "", strsplit(filestr, "\\} \\{")[[1]])
    else
        files <- strsplit(filestr, " ")[[1]]

    if(!is.na(language))
        language <- paste("\"", language, "\"", sep="")

    if(encoding == "")
        encoding <- "auto"

    if(length(files) == 1) {
        doItAndPrint(sprintf('corpus <- Corpus(AlcesteSource("%s", "%s"), readerControl=list(language=%s))',
                             files[1], encoding, language))

		if(!exists("corpus") || length(corpus) == 0) {
		    .Message(.gettext("Reading the specified file failed. Are you sure this file is in the correct format?"),
		             type="error")

		    return(FALSE)
		}
    }
    else {
        doItAndPrint(sprintf('corpusList <- vector(%s, mode="list")', length(files)))
        for(i in seq(length(files))) {
            doItAndPrint(sprintf('corpusList[[%s]] <- Corpus(AlcesteSource("%s"), readerControl=list(language=%s))',
                                 i, files[i], language))

			if(length(corpusList[[i]]) == 0) {
				.Message(.gettext("Reading the specified file failed. Are you sure this file is in the correct format?"),
				         type="error")

				return(FALSE)
			}
        }

        doItAndPrint("corpus <- do.call(c, c(corpusList, list(recursive=TRUE)))")
        doItAndPrint("rm(corpusList)")
        gc()
    }

    # We rely on names/IDs later e.g. in showCorpusCa() because we cannot use indexes when documents are skipped
    # Duplicated IDs can happen since they can be specified in the Alceste file, but also set automatically
    # when missing
    doItAndPrint("names(corpus) <- make.unique(names(corpus))")

    doItAndPrint("corpusVars <- extractMetadata(corpus, date=FALSE)")

    list(source=sprintf(.ngettext(length(files), "Alceste file %s", "Alceste files %s"),
                        paste(files, collapse=", ")))
}

# Choose a Twitter hashtag to search for messages
importCorpusFromTwitter <- function(language=NA) {
    if(!.checkAndInstall("twitteR",
                         .gettext("The twitteR package is needed to import corpora from Twitter.\nDo you want to install it?")))
        return(FALSE)

    if(!is.na(language))
        language <- paste("\"", language, "\"", sep="")

    initializeDialog(title=.gettext("Import Corpus From Twitter"))

    tclReqURL <- tclVar("https://api.twitter.com/oauth/request_token")
    entryReqURL <- ttkentry(top, width=37, textvariable=tclReqURL)

    tclAuthURL <- tclVar("https://api.twitter.com/oauth/authorize")
    entryAuthURL <- ttkentry(top, width=37, textvariable=tclAuthURL)

    tclAccessURL <- tclVar("https://api.twitter.com/oauth/access_token")
    entryAccessURL <- ttkentry(top, width=37, textvariable=tclAccessURL)

    tclConsumerKey <- tclVar("")
    entryConsumerKey <- ttkentry(top, width=37, textvariable=tclConsumerKey)

    tclConsumerSecret <- tclVar("")
    entryConsumerSecret <- ttkentry(top, width=37, textvariable=tclConsumerSecret)

    # TRANSLATORS: replace 'en' with your language's ISO 639 two-letter code
    tclText <- tclVar("")
    entryText <- ttkentry(top, width=37, textvariable=tclText)

    tclNMess <- tclVar(100)
    tclNSlider <- tkscale(top, from=1, to=1500,
                          showvalue=TRUE, variable=tclNMess,
                          resolution=1, orient="horizontal")

    checkBoxes(frame="optionsFrame",
               boxes=c("removeNames", "removeHashtags", "exclRetweets"),
               initialValues=c(1, 1, 0),
               labels=c(.gettext("Remove user names"), .gettext("Remove hashtags"),
                        .gettext("Exclude retweets")),
               title=.gettext("Options:"))

    result <- tclVar()

    onOK <- function() {
        reqURL <- tclvalue(tclReqURL)
        authURL <- tclvalue(tclAuthURL)
        accessURL <- tclvalue(tclAccessURL)
        consumerKey <- tclvalue(tclConsumerKey)
        consumerSecret <- tclvalue(tclConsumerSecret)
        text <- tclvalue(tclText)
        nmess <- tclvalue(tclNMess)
        exclRetweets <- tclvalue(exclRetweetsVariable) == 1

        if(reqURL == "" || authURL == "" || accessURL == "" ||
           consumerKey == "" || consumerSecret == "") {
            .Message(.gettext("Please enter valid authentication settings."), type="error", parent=top)
            return(FALSE)
        }
        if(text == "") {
            .Message(.gettext("Please enter valid text to search for."), type="error", parent=top)
            return(FALSE)
        }

        closeDialog()

        setBusyCursor()
        on.exit(setIdleCursor())

        # In case something goes wrong
        tclvalue(result) <- "error"

        doItAndPrint("library(twitteR)")

        doItAndPrint(sprintf('twitCred <- OAuthFactory$new(consumerKey="%s", consumerSecret="%s", requestURL="%s", accessURL="%s", authURL="%s")', consumerKey, consumerSecret, reqURL, accessURL, authURL))
        doItAndPrint("twitCred$handshake()")

        if(!isTRUE(twitCred$handshakeComplete)) {
            .Message(.gettext("TwitteR authentication failed. Please check the entered credentials or PIN code."),
                     type="error", parent=top)
            return(FALSE)
        }

        doItAndPrint("registerTwitterOAuth(twitCred)")

        doItAndPrint(sprintf('messages <- searchTwitter("%s", %s, %s)', text, nmess, language))

        if(length(messages) == 0) {
            .Message(sprintf(.gettext("No recent tweets match the specified search criteria in the chosen language (%s)."), language),
                     type="error", parent=top)
            return(FALSE)
        }

        doItAndPrint("corpusDataset <- twListToDF(messages)")

        if(length(unique(strftime(corpusDataset$created, "%y-%m-%d"))) == 1)
            fmt <- "%H:%M"
        else if(length(unique(strftime(corpusDataset$created, "%y-%m"))) == 1)
            fmt <- "%d %H:%M"
        else if(length(unique(strftime(corpusDataset$created, "%y"))) == 1)
            ftm <- "%m-%d %H:%M"
        else
            fmt <- "%y-%m-%d %H:%M"

        doItAndPrint(sprintf('rownames(corpusDataset) <- make.unique(paste(abbreviate(corpusDataset$screenName, 10), strftime(corpusDataset$created, "%s")))', fmt))

        doItAndPrint(sprintf('corpus <- Corpus(DataframeSource(corpusDataset[1]), readerControl=list(language=%s))',
                             language))
        doItAndPrint("rm(messages)")

        if(!exists("corpus") || length(corpus) == 0) {
            .Message(.gettext("Retrieving messages from Twitter failed."),
                     type="error", parent=top)
            return(FALSE)
        }

        doItAndPrint('corpusVars <- corpusDataset[c("screenName", "created", "truncated", "statusSource")]')
        doItAndPrint("rm(corpusDataset)")
        doItAndPrint(sprintf('colnames(corpusVars) <- c("%s", "%s", "%s", "%s")',
                             # Keep in sync with extractMetadata()
                             .gettext("Author"), .gettext("Time"), .gettext("Truncated"), .gettext("StatusSource")))

        doItAndPrint(sprintf('corpusVars[["%s"]] <- grepl("\\\\bRT\\\\b", corpus)', .gettext("Retweet")))

        if(exclRetweets) {
            doItAndPrint(sprintf('corpus <- corpus[!corpusVars[["%s"]]]', .gettext("Retweet")))
            doItAndPrint(sprintf('corpusVars <- subset(corpusVars, !%s)', .gettext("Retweet")))
        }

        tclvalue(result) <- "success"

        return(FALSE)
    }

    onCancel <- function() {
        if (GrabFocus()) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(CommanderWindow())
        tclvalue(result) <- "cancel"
    }

    OKCancelHelp(helpSubject="importCorpusDlg")
    tkgrid(labelRcmdr(top, text=.gettext("Note: Twitter requires you to register a custom application and fill in\nthe details below. See vignette(\"twitteR\") and https://dev.twitter.com/apps/new/.\nYou will need to switch manually to the R console and copy the PIN\ncode you get from the URL printed there.")), sticky="w", pady=6, columnspan=2)
    tkgrid(labelRcmdr(top, text=.gettext("Request token URL:")),
           entryReqURL, sticky="w", pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Authorize URL:")),
           entryAuthURL, sticky="w", pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Access token URL:")),
           entryAccessURL, sticky="w", pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Consumer key:")),
           entryConsumerKey, sticky="w", pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Consumer secret:")),
           entryConsumerSecret, sticky="w", pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Text to search for:")),
           entryText, sticky="w", pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Maximum number of tweets to download:")),
           tclNSlider, sticky="w", pady=6)
    tkgrid(optionsFrame, sticky="w", pady=6, columnspan=2)
    tkgrid(buttonsFrame, columnspan=2, sticky="ew", pady=6)
    dialogSuffix(focus=entryText, force.wait=TRUE)

    if(tclvalue(result) == "success")
        return(list(source=sprintf(.gettext("Twitter search for %s"), tclvalue(tclText)),
                    removeNames=tclvalue(removeNamesVariable) == 1,
                    removeHashtags=tclvalue(removeHashtagsVariable) == 1))
    else
        return(FALSE)
}

# Adapted version of tm's makeChunks() remembering which chunk comes from which document,
# preserving corpus meta-data, and skipping empty chunks.
# Copyright Ingo Feinerer, Licence: GPL (>= 2).
# http://tm.r-forge.r-project.org/
splitTexts <- function (corpus, chunksize, preserveMetadata=TRUE) 
{
    chunks <- list(length(corpus))
    origins <- list(length(corpus))

    for (k in seq_along(corpus)) {
        chunks_k <- tapply(content(corpus[[k]]),
                           rep(seq(1, length(content(corpus[[k]]))),
                               each=chunksize, length.out=length(content(corpus[[k]]))), c)

        # Skeep empty chunks
        keep <- nchar(gsub("[\n[:space:][:punct:]]+", "", sapply(chunks_k, paste, collapse=""))) > 0

        chunks[[k]] <- chunks_k[keep]
        origins[[k]] <- rep(k, sum(keep))
    }

    # Merge only the per-document lists of chunks at the end to reduce the number of copies
    chunks <- do.call(c, chunks)
    origins <- do.call(c, origins)

    newCorpus <- Corpus(VectorSource(chunks))

    names1 <- names(corpus)
    names2 <- make.unique(names1[origins])

    # Copy meta data from old documents
    if(preserveMetadata) {
        newCorpus$dmeta <- meta(corpus)[origins,, drop=FALSE]

        for(i in seq_along(corpus)) {
            attrs <- meta(corpus[[i]])

            for(j in which(origins == i)) {
                doc <- newCorpus[[j]]
                doc$meta <- attrs
                meta(doc, "id") <- names2[j]
                meta(doc, "document") <- names1[i]
                newCorpus[[j]] <- doc
            }
        }
    }

    meta(newCorpus, .gettext("Doc ID")) <- names1[origins]
    meta(newCorpus, .gettext("Doc N")) <- origins
    names(newCorpus) <- names2

    newCorpus
}
