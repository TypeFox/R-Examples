printLanguages <- function(){
    # Convenience user-facing function to print language
    # codes for HP API.
    #    
    # Args:
    #   NONE
    #    
    # Returns:
    #   NOTHING. Prints list of supported languages
    #   and codes. 
    language.codes <- list(
        German   = "de-DE",
        USA      = "en-US",	
        GB       = "en-GB",	
        Spanish  = "es-ES",	
        French   = "fr-FR",	
        Italian  = "it-IT",	
        Mandarin = "zh-CN"	
        )
    cat('LANGUAGE                 CODE
-----------------------------
')
    for(i in 1:length(language.codes)){
        code <- language.codes[i]
        cat(paste(gsub('_', ' ', names(code)), paste(rep(" ", 22 - nchar(names(code))), collapse=''), code))
        cat('\n')
    }
}
