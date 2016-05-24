translate <-
function(dataset = NULL, content.field = NULL, content.vec = NULL,
                      google.api.key = NULL, microsoft.client.id = NULL, microsoft.client.secret = NULL,
                      source.lang = NULL, target.lang = NULL){

    # Do some sanity checking
    translator <- validateInput(dataset, content.field, content.vec, google.api.key,
                                microsoft.client.id, microsoft.client.secret,
                                source.lang, target.lang)
    
    # Get translation vector
    if(!(is.null(dataset))){
        to.translate <- dataset[[content.field]]
    }
    if(!(is.null(content.vec))){
        to.translate <- content.vec
    }

    checkLang(to.translate, source.lang, translator)
    
    # Do the translation
    if(translator == 'Google'){
        translated <- unname(
            unlist(
                mclapply(to.translate, function(x) googleTranslate(x, google.api.key, source.lang, target.lang))
                )
            )
    }

    if(translator == 'Microsoft'){ 
        ptm <- proc.time()
        access.token <- getAccessToken(microsoft.client.id, microsoft.client.secret)
        translated <- c()
        for(doc in to.translate){
            translated <- c(translated, microsoftTranslate(doc, access.token, source.lang, target.lang))
            if((proc.time() - ptm)[3] > 540){
                ptm <- proc.time()
                access.token <- getAccessToken(microsoft.client.id, microsoft.client.secret)
            }
        }
    }

    # Figure out what we should return
    if(!(is.null(content.vec))){
        return(translated)
    }
    if(!(is.null(dataset) & is.null(content.field))){
        dataset$translatedContent <- translated
        return(dataset)
    }
}
