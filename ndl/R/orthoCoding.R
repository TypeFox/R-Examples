orthoCoding = function (strings=c("hel.lo","wor.ld"), grams = c(2), tokenized = F, sepToken = '.') {
    
    if (length(grams) < 1) { stop("This function requires a non-zero length vector of n-gram sizes for the argument 'grams'.")}

    if (! is.numeric(grams)) { stop("This function requires a vector of one or more numbers of n-gram sizes for the argument 'grams'.")}
    
    ngram.fnc = function(s, n) {
        if (n == 1) { # remove hash cues for unigrams
            s = sub(paste('#',sepToken,sep=''),'',s)
            s = sub(paste(sepToken,'#',sep=''),'',s)
        }
        tokens = unlist(strsplit(s, sepToken, fixed=T))
        len = length(tokens)
        ng = NULL
        for (i in 1:(len - n + 1)) {
            ng = c(ng, paste(tokens[i:(i+n-1)],collapse=''))
        }
        return(paste(ng, collapse = "_"))
    }

    if (!tokenized) sepToken = ''
    
    letters = sapply(strings, FUN = function(s) paste('#',sepToken, s,
sepToken,'#', sep = ""))

    for (i in 1:length(grams)) {
        cuesi = unlist(lapply(letters, FUN = ngram.fnc, grams[i]))
        if (exists("mycues") == 0) {
            mycues = cuesi
        } else {
            mycues = paste(mycues, cuesi, sep = "_")
        }
    }
    return(mycues)
}

## Old version... now depreciated.
## orthoCoding <-
## function (words = c("hello", "world"), maxn = 2, inclusive=FALSE) 
## {
##     ngram.fnc = function(s, n) {
##         len = nchar(s)
##         ng = NULL
##         for (i in 1:(len - n + 1)) {
##             ng = c(ng, substr(s, i, i + n - 1))
##         }
##         return(paste(ng, collapse = "_"))
##     }
##     letters = strsplit(words, "")
##     grams = unlist(lapply(letters, FUN = paste, collapse = "_"))
##     letters = sapply(words, FUN = function(s) paste("#", s, "#", sep = ""))
##     if (maxn == 1) {
##         return(grams)
##     }
##     else {
##       if (inclusive) {
##         for (i in 2:maxn) {
##           gramsi = unlist(lapply(letters, FUN = ngram.fnc, 
##             i))
##           grams = paste(grams, gramsi, sep = "_")
##         }
##       } else {
##         grams = unlist(lapply(letters, FUN = ngram.fnc, maxn))
##         return(grams)        
##       }
##     }
##     return(grams)
##   }
