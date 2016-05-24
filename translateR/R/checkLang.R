checkLang <-
function(to.translate, source.lang, translator){
    languages <- languageCodes()
    combined.vec <- paste(to.translate, collapse = ' ')
    guessed.lang <- toupper(textcat(combined.vec))
    indicated.lang <- toupper(names(languages[[translator]][unlist(languages[[translator]]) == source.lang]))

    if(!(guessed.lang %in% indicated.lang)){
        msg <- paste("\nThe content appears to be in ", guessed.lang, ". However, the language code you provided suggests that the text is in ", indicated.lang, ". If you entered the wrong language code, stop the process. Otherwise, translateR will treat the text as ", indicated.lang, '.', sep = '')

        msg = strwrap(msg, width = 0.9 * getOption("width"))
        msg <- paste(msg, collapse="\n")       
        warning(msg, call. = FALSE)
    }
}
