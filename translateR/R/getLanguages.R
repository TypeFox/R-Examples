getLanguages <-
function(translator){
    languages <- languageCodes()
    
    cat('LANGUAGE                CODE
----------------------------
')
    for(i in 1:length(languages[[translator]])){
        code <- languages[[translator]][i]
        cat(paste(gsub('_', ' ', names(code)), paste(rep(" ", 22 - nchar(names(code))), collapse=''), code))
        cat('\n')
    }
}
