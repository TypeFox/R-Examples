#' @demoTitle dallas-text
#' 
#' Demo text tf, tf-idf and parser functions
#'
#' To install and use Dallas demo dataset in Aster:
#'
#' 1. download DallasOpenData.zip from
#'   https://bitbucket.org/grigory/toaster/downloads/DallasOpenData.zip
#' 2. run script to create data set in Aster
#'   sh load_dallas_data.sh -d mydbname -U username -w mypassword 
#' 3. create Aster ODBC DSN on your desktop
#'   see https://bitbucket.org/grigory/toaster/wiki/Home#markdown-header-odbc-driver-and-dns

library(toaster)
library(tm)
library(RColorBrewer)
library(reshape2)

# utility input function
readlineDef <- function(prompt, default) {
  if (!is.null(prompt))
    prompt = paste0(prompt, "[", default, "]: ")
  else 
    prompt = paste0(prompt, ": ")
  
  result = readline(prompt)
  if (result == "") 
    return (default)
  else
    return (result)
}

# utility connection function
connectWithDSNToAster <- function(dsn=NULL) {
  dsn = readlineDef("Enter Aster ODBC DSN: ", dsn)
  
  tryCatch(close(conn), error=function(err) {NULL})
  
  conn = tryCatch({
    conn = odbcConnect(dsn)
    odbcGetInfo(conn)
    return (conn)
  }, error=function(err) {
    stop(paste("Can't connect to Aster - check DSN '", dsn, "'"))
  })
}

## connect to Aster first
conn = connectWithDSNToAster()

# term frequency on narratives
# divided into 4 time groups: 0-3
tfIdfDallasByTime = computeTfIdf(channel=conn, tableName="public.dallaspoliceall", 
                                 docId="(extract('hour' from offensestarttime)::int/6)%4", 
                    textColumns=c('offensedescription','offensenarrative'),
                    parser=nGram(1, punctuation="[\\\\[\\\\]%/`~#^&*()''\":;><@0-9-]+"),
                    test=FALSE)
dim(inspect(tfIdfDallasByTime))

findFreqTerms(tfIdfDallasByTime, 0.00005)

findAssocs(tfIdfDallasByTime, "INVESTIGATIONS", 0.99)

# top terms (by tf) for each document
tfIdfDallasByTimeDF <- cbind(term=dimnames(tfIdfDallasByTime)[[1]], as.data.frame(inspect(tfIdfDallasByTime)), stringsAsFactors=FALSE)
tfIdfDallasByTimeDFMelt = melt(tfIdfDallasByTimeDF, id.vars='term', measure.vars=as.character(0:3), 
                          variable.name='docid', value.name='tf')
tfOrdered <- tfIdfDallasByTimeDFMelt[order(-tfIdfDallasByTimeDFMelt$tf), ]
by(tfOrdered, tfOrdered["docid"], head, n=20)

## define stop words for post-processing
stopwords = unique(c(letters, "a", "about", "above", "above", "across", "after", "afterwards", "again", "against", "all", "almost", "alone", "along", "already", "also","although","always","am","among", "amongst", "amoungst", "amount",  "an", "and", "another", "any","anyhow","anyone","anything","anyway", "anywhere", "are", "around", "as",  "at", "back","be","became", "because","become","becomes", "becoming", "been", "before", "beforehand", "behind", "being", "below", "beside", "besides", "between", "beyond", "bill", "both", "bottom","but", "by", "call", "can", "cannot", "cant", "co", "con", "could", "couldnt", "cry", "de", "describe", "detail", "do", "done", "down", "due", "during", "each", "eg", "eight", "either", "eleven","else", "elsewhere", "empty", "enough", "etc", "even", "ever", "every", "everyone", "everything", "everywhere", "except", "few", "fifteen", "fify", "fill", "find", "fire", "first", "five", "for", "former", "formerly", "forty", "found", "four", "from", "front", "full", "further", "get", "give", "go", "had", "has", "hasnt", "have", "he", "hence", "her", "here", "hereafter", "hereby", "herein", "hereupon", "hers", "herself", "him", "himself", "his", "how", "however", "hundred", "ie", "if", "in", "inc", "indeed", "interest", "into", "is", "it", "its", "itself", "keep", "last", "latter", "latterly", "least", "less", "ltd", "made", "many", "may", "me", "meanwhile", "might", "mill", "mine", "more", "moreover", "most", "mostly", "move", "much", "must", "my", "myself", "name", "namely", "neither", "never", "nevertheless", "next", "nine", "no", "nobody", "none", "noone", "nor", "not", "nothing", "now", "nowhere", "of", "off", "often", "on", "once", "one", "only", "onto", "or", "other", "others", "otherwise", "our", "ours", "ourselves", "out", "over", "own","part", "per", "perhaps", "please", "put", "rather", "re", "same", "see", "seem", "seemed", "seeming", "seems", "serious", "several", "she", "should", "show", "side", "since", "sincere", "six", "sixty", "so", "some", "somehow", "someone", "something", "sometime", "sometimes", "somewhere", "still", "such", "system", "take", "ten", "than", "that", "the", "their", "them", "themselves", "then", "thence", "there", "thereafter", "thereby", "therefore", "therein", "thereupon", "these", "they", "thickv", "thin", "third", "this", "those", "though", "three", "through", "throughout", "thru", "thus", "to", "together", "too", "top", "toward", "towards", "twelve", "twenty", "two", "un", "under", "until", "up", "upon", "us", "very", "via", "was", "we", "well", "were", "what", "whatever", "when", "whence", "whenever", "where", "whereafter", "whereas", "whereby", "wherein", "whereupon", "wherever", "whether", "which", "while", "whither", "who", "whoever", "whole", "whom", "whose", "why", "will", "with", "within", "without", "would", "yet", "you", "your", "yours", "yourself", "yourselves", "the"))

## use regex to remove numbers using stop word list
stopwords = c(stopwords, "\\d+")

tfIdfDallasByRace = computeTfIdf(channel=conn, tableName="public.dallaspoliceall", 
                docId="offenserace", 
                textColumns=c('offensedescription','offensenarrative'),
                parser=nGram(2, punctuation="[\\\\[\\\\]%/`~#^&*()-]+"),
                where="offenserace in ('B','L','A','W','C')", stopwords=stopwords)
words = with(tfIdfDallasByRace$rs, tfIdfDallasByRace$rs[docid=='W',])
pal = rev(brewer.pal(8, "Set1"))[c(-3,-1)]
createWordcloud(words$term, words$tf, maxWords=20, scale=c(4,1), palette=pal)

tfidf1 = computeTfIdf(channel=conn, tableName="public.dallaspoliceall", docId="offenserace", 
                      textColumns=c('offensedescription','offensenarrative'),
                      parser=nGram(2), where="offenserace in ('B','W','L','A','C')",
                      stopwords=stopwords)
words=with(tfidf1$rs, tfidf1$rs[docid=='A',])
createWordcloud(words$term, words$tf_idf, maxWords=25, scale=c(4, 0.5), palette=pal, 
                file='offense_asian.png', width=700, height=700)

tfidf1 = computeTfIdf(channel=conn, tableName="texts.containertrailerplanpaths", docId="orig", 
                      textColumns=c("f5466all", "commentsall"),
                      where="f5466all not like '%null%' and commentsall not like '%null%'",
                      parser=nGram(2:3, punctuation="[\\\\[\\\\]%/`~#^&*()-]+"))

sort(table(tfidf1$docid))
words = tfidf1[tfidf1$docid == '913',]
createWordcloud(words$term, words$tf_idf, maxWords=25, scale=c(4, 0.5), title="Top 25 2-grams for 913")



