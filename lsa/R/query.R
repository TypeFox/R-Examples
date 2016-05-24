### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### query.R
### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### dependencies: library("RStem")
### 
### 2014-03-28: replaced stemming from Snowball functions to SnowballC functions
### 2009-08-17: replaced wordStem with SnowballStemmer
### 2005-11-12: bugfix for regexp
### 2005-11-08: pseudo_dtm.R renamed to query.R
### 2005-11-08: removed pseudo_svd (now in lsa.R, renamed
###             to foldinLSAspace()
### 2005-11-08: removed pseudo_docs function (can be done
###             with textmatrix, now!)
### 2005-08-25: added "\\[|\\]|\\{|\\}" to gsub
### 

query <- function( qtext, termlist, stemming=FALSE, language="german" ) {
    
    # qtext: string with the query words, whitespace separated
    # termlist: list of allowed terms
    # dtm: original doc-term-matrix (no weighting applied!)
    
    dtm = NULL
    
    q = strsplit( gsub('[[:space:]]|[[:punct:]]+', ' ', tolower(qtext) ), " ")[[1]]
    vec = vector( mode="numeric", length(termlist) )
    for ( word in q ) {
        if (stemming) word = wordStem(word, language)
        if (word != "") {
            vec[ match(word,termlist) ] = vec[ match(word,termlist) ] + 1
        }
    }
    
    dtm = as.matrix(vec)
    colnames(dtm) = toupper(qtext)
    rownames(dtm) = termlist
    
    environment(dtm) = new.env()
    class(dtm) = "textmatrix"
    
    return ( dtm )
    
}

