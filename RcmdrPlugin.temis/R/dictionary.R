termsDictionaryAlpha <- function() {
    doItAndPrint('dict <- termsDictionary(dtm)')

    setLastTable("dict", .gettext("Terms dictionary sorted alphabetically"))

    doItAndPrint("dict")

    activateMenus()
}

termsDictionaryOcc <- function() {
    doItAndPrint('dict <- termsDictionary(dtm, "occurrences")')

    setLastTable("dict", .gettext("Terms dictionary sorted by number of occurrences"))

    doItAndPrint("dict")

    activateMenus()
}

termsDictionary <- function(dtm, order=c("alphabetic", "occurrences")) {
    order <- match.arg(order)

    setBusyCursor()
    on.exit(setIdleCursor())

    lang <- attr(dtm, "language")
    processing <- attr(dtm, "processing")
    dict <- attr(dtm, "dictionary")
    stopword <- rownames(dict) %in% stopwords(lang)

    if(processing["stemming"] || processing["customStemming"]) {
        dict <- cbind(dict[-3],
                      col_sums(dtm)[dict[[2]]],
                      dict[3],
                      # Some words can be removed as stopwords, but be present because another
                      # word that has been kept is identical in its stemmed from
                      ifelse(!dict[[2]] %in% Terms(dtm) |
                             ((!processing["customStemming"] & processing["stopwords"] & stopword) |
                              (processing["customStemming"] & dict[[2]] == "")),
                             .gettext("Removed"), ""))

        colnames(dict)[3:5] <- c(.gettext("Stemmed occ."), .gettext("Stopword"), .gettext("Removed"))

        if(order == "occurrences")
            dict[order(dict[[3]], decreasing=TRUE),]
        else
            # When editing stemming manually, words to exclude are set to ""
            dict[order(ifelse(dict[[2]] == "", rownames(dict), dict[[2]])),]
    }
    else {
        dict <- cbind(dict,
                      # Some words can be removed as stopwords, but be present because another
                      # word that has been kept is identical in its stemmed from
                      ifelse(!rownames(dict) %in% Terms(dtm) | (stopword & processing["stopwords"]),
                             .gettext("Removed"), ""))
        colnames(dict)[2:3] <- c(.gettext("Stopword"), .gettext("Removed"))

        if(order == "occurrences")
            dict[order(dict[[1]], decreasing=TRUE),]
        else
            dict

    }
}
