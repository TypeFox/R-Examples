#' Compute term frequencies on a corpus.
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param tableName Aster table name
#' @param docId vector with one or more column names comprising unique document id. 
#'   Values are concatenated with \code{idSep}. Database NULLs are replaced with
#'   \code{idNull} string.
#' @param textColumns one or more names of columns with text. Multiple coumn are
#'   concatenated into single text field first.
#' @param parser type of parser to use on text. For example, \code{ngram(2)} parser
#'   generates 2-grams (ngrams of length 2), \code{token(2)} parser generates 2-word 
#'   combinations of terms within documents.
#' @param weighting term frequency formula to compute the tf value. One of following: 
#'   \code{'raw'}, \code{'bool'}, \code{'binary'}, \code{'log'}, \code{'augment'}, and
#'   \code{'normal'} (default). 
#' @param top specifies threshold to cut off terms ranked below \code{top} value. If value
#'   is greater than 0 then included top ranking terms only, otherwise all terms returned 
#'   (also see paramter \code{rankFunction}). Terms are always ordered by their term frequency (tf)
#'   within each document. Filtered out terms have their rank ariphmetically greater than 
#'   threshold \code{top} (see details): term is more important the smaller value of its rank.
#' @param rankFunction one of \code{rownumber, rank, denserank, percentrank}. Rank computed and
#'   returned for each term within each document. function determines which SQL window function computes 
#'   term rank value (default \code{rank} corresponds to SQL \code{RANK()} window function). 
#'   When threshold \code{top} is greater than 0 ranking function used to limit number of 
#'   terms returned (see details).
#' @param idSep separator when concatenating 2 or more document id columns (see \code{docId}).
#' @param idNull string to replace NULL value in document id columns.
#' @param where specifies criteria to satisfy by the table rows before applying
#'   computation. The criteria are expressed in the form of SQL predicates (inside 
#'   \code{WHERE} clause).
#' @param stopwords character vector with stop words. Removing stop words takes place in R after 
#'   results are computed and returned from Aster.
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \link{RODBC} 
#'   functions \link{sqlQuery} and \link{sqlSave}). 
#'   
#' @details
#' By default function computes and returns all terms. When large number of terms is expected then
#' use parameters \code{top} to limit number of terms returned by 
#' filtering top ranked terms for each document. Thus if set \code{top=1000} and there
#' is 100 documents then at least 100,000 terms (rows) will be returned. Result size could 
#' exceed this number when other than \code{rownumber} \code{rankFunction} used:
#' \itemize{
#'     \item \emph{\code{rownumber}} applies a sequential row number, starting at 1, to each term in a document.
#'       The tie-breaker behavior is as follows: Rows that compare as equal in the sort order will be
#'       sorted arbitrarily within the scope of the tie, and all terms will be given unique row numbers.
#'     \item \emph{\code{rank}} function assigns the current row-count number as the terms's rank, provided the 
#'       term does not sort as equal (tie) with another term. The tie-breaker behavior is as follows: 
#'       terms that compare as equal in the sort order are sorted arbitrarily within the scope of the tie, 
#'       and the sorted-as-equal terms get the same rank number.
#'     \item \emph{\code{denserank}} behaves like the \code{rank} function, except that it never places 
#'       gaps in the rank sequence. The tie-breaker behavior is the same as that of RANK(), in that 
#'       the sorted-as-equal terms receive the same rank. With \code{denserank}, however, the next term after 
#'       the set of equally ranked terms gets a rank 1 higher than preceding tied terms.
#'     \item \emph{\code{percentrank}} assigns a relative rank to each term, using the formula: 
#'       \code{(rank - 1) / (total rows - 1)}. The tie-breaker behavior is as follows: Terms that compare 
#'       as equal are sorted arbitrarily within the scope of the tie, and the sorted-as-equal rows 
#'       get the same percent rank number.
#' }
#' The ordering of the rows is always by their tf value within each document.
#'  
#' @seealso \code{computeTfIdf}, \code{\link{nGram}}, \code{\link{token}}
#' @export 
#' @examples
#' if(interactive()){
#' # initialize connection to Dallas database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' # compute term-document-matrix of all 2-word Ngrams of Dallas police open crime reports
#' tdm1 = computeTf(channel=conn, tableName="public.dallaspoliceall", docId="offensestatus",
#'                  textColumns=c("offensedescription", "offensenarrative"),
#'                  parser=nGram(2),
#'                  where="offensestatus NOT IN ('System.Xml.XmlElement', 'C')")
#'
#' # compute term-document-matrix of all 2-word combinations of Dallas police crime reports
#' # by time of day (4 documents corresponding to 4 parts of day)
#' tdm2 = computeTf(channel=conn, tableName="public.dallaspoliceall",
#'                  docId="(extract('hour' from offensestarttime)/6)::int%4",
#'                  textColumns=c("offensedescription", "offensenarrative"),
#'                  parser=token(2, punctuation="[-.,?\\!:;~()]+", stopWords=TRUE),
#'                  where="offensenarrative IS NOT NULL")
#'
#' # include only top 100 ranked 2-word ngrams for each offense status
#' # into resulting term-document-matrix using dense rank function
#' tdm3 = computeTf(channel=NULL, tableName="public.dallaspoliceall", docId="offensestatus",
#'                  textColumns=c("offensedescription", "offensenarrative"),
#'                  parser=nGram(2), top=100, rankFunction="denserank",
#'                  where="offensestatus NOT IN ('System.Xml.XmlElement', 'C')")
#' 
#' }
#' 
computeTf <- function(channel, tableName, docId, textColumns, parser,
                      weighting = "normal",
                      top = NULL, rankFunction = "rank",
                      where = NULL, idSep = '-', idNull = '(null)',
                      stopwords = NULL, test = FALSE) {
  
  weighting = match.arg(weighting, c('raw','bool','binary','log','augment','normal'))
  rankFunction = match.arg(rankFunction, c('rank', 'rownumber', 'row', 'denserank', 'percentrank'))
  tfFormula = switch(tolower(weighting),
                    normal="normal",
                    raw="normal",
                    bool="bool",
                    binary="bool",
                    log="log",
                    augment="augment"
  )
  
  if (missing(tableName)) {
    stop("Table name must be specified.")
  }
  
  if (missing(docId)) {
    stop("Doc id must be specified.")
  }
  
  if (missing(textColumns) || length(textColumns)==0) {
    stop("Text columns must be specified.")
  }
  
  isValidConnection(channel, test)
  
  windowFunction = getWindowFunction(rankFunction)
  
  where_clause = makeWhereClause(where)
  
  derivedDocId = makeDocumentId(docId, idSep, idNull)
  
  textSelectSQL = parseTextSQL(parser, tableName, derivedDocId, textColumns, where)
  
  sql = paste0(
    "SELECT * FROM 
       (SELECT *, ", windowFunction, " OVER (PARTITION BY docid ORDER BY tf DESC) rank FROM TF(
         ON (SELECT docid, term FROM ( ", textSelectSQL, " ) t ) PARTITION BY docid
         FORMULA('", tfFormula, "')
       )) t2", makeRankFilter(top)
    )
  
  if (test) 
    return (sql)
  else {
    rs = toaSqlQuery(channel, sql, stringsAsFactors = FALSE, na.strings="")
  }
  
  rs = removeStopWords(rs, stopwords)
  
  x = makeSimpleTripletMatrix(rs, ifelse(weighting == 'raw', 'count', 'tf'), 'tf')
  return(x)
}

#' Compute Term Frequency - Inverse Document Frequency on a corpus.
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param tableName Aster table name
#' @param docId vector with one or more column names comprising unique document id. 
#'   Values are concatenated with \code{idSep}. Database NULLs are replaced with
#'   \code{idNull} string.
#' @param textColumns one or more names of columns with text. Multiple coumn are
#'   concatenated into single text field first.
#' @param parser type of parser to use on text. For example, \code{ngram(2)} parser
#'   generates 2-grams (ngrams of length 2), \code{token(2)} parser generates 2-word 
#'   combinations of terms within documents.
#' @param top specifies threshold to cut off terms ranked below \code{top} value. If value
#'   is greater than 0 then included top ranking terms only, otherwise all terms returned 
#'   (also see paramter \code{rankFunction}). Terms are always ordered by their term frequency -
#'   inverse document frequency (tf-idf) within each document. Filtered out terms have their 
#'   rank ariphmetically greater than threshold \code{top} (see details): term is more 
#'   important the smaller value of its rank.
#' @param rankFunction one of \code{rownumber, rank, denserank, percentrank}. Rank computed and
#'   returned for each term within each document. function determines which SQL window function computes 
#'   term rank value (default \code{rank} corresponds to SQL \code{RANK()} window function). 
#'   When threshold \code{top} is greater than 0 ranking function used to limit number of 
#'   terms returned (see details).
#' @param idSep separator when concatenating 2 or more document id columns (see \code{docId}).
#' @param idNull string to replace NULL value in document id columns.
#' @param adjustDocumentCount logical: if TRUE then number of documents 2 will be increased by 1.
#' @param where specifies criteria to satisfy by the table rows before applying
#'   computation. The criteria are expressed in the form of SQL predicates (inside 
#'   \code{WHERE} clause).
#' @param stopwords character vector with stop words. Removing stop words takes place in R after 
#'   results are computed and returned from Aster.
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \link{RODBC} 
#'   functions \link{sqlQuery} and \link{sqlSave}). 
#'   
#' @details
#' By default function computes and returns all terms. When large number of terms is expected then
#' use parameters \code{top} to limit number of terms returned by 
#' filtering top ranked terms for each document. Thus if set \code{top=1000} and there
#' is 100 documents then at least 100,000 terms (rows) will be returned. Result size could 
#' exceed this number when other than \code{rownumber} \code{rankFunction} used:
#' \itemize{
#'     \item \emph{\code{rownumber}} applies a sequential row number, starting at 1, to each term in a document.
#'       The tie-breaker behavior is as follows: Rows that compare as equal in the sort order will be
#'       sorted arbitrarily within the scope of the tie, and all terms will be given unique row numbers.
#'     \item \emph{\code{rank}} function assigns the current row-count number as the terms's rank, provided the 
#'       term does not sort as equal (tie) with another term. The tie-breaker behavior is as follows: 
#'       terms that compare as equal in the sort order are sorted arbitrarily within the scope of the tie, 
#'       and the sorted-as-equal terms get the same rank number.
#'     \item \emph{\code{denserank}} behaves like the \code{rank} function, except that it never places 
#'       gaps in the rank sequence. The tie-breaker behavior is the same as that of RANK(), in that 
#'       the sorted-as-equal terms receive the same rank. With \code{denserank}, however, the next term after 
#'       the set of equally ranked terms gets a rank 1 higher than preceding tied terms.
#'     \item \emph{\code{percentrank}} assigns a relative rank to each term, using the formula: 
#'       \code{(rank - 1) / (total rows - 1)}. The tie-breaker behavior is as follows: Terms that compare 
#'       as equal are sorted arbitrarily within the scope of the tie, and the sorted-as-equal rows 
#'       get the same percent rank number.
#' }
#' The ordering of the rows is always by their tf-idf value within each document.   
#'   
#' @seealso \code{computeTf}, \code{\link{nGram}}, \code{\link{token}}
#' @export 
#' @examples
#' if(interactive()){
#' # initialize connection to Dallas database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' # compute term-document-matrix of all 2-word Ngrams of Dallas police crime reports
#' # for each 4-digit zip
#' tdm1 = computeTfIdf(channel=conn, tableName="public.dallaspoliceall", 
#'                     docId="substr(offensezip, 1, 4)", 
#'                     textColumns=c("offensedescription", "offensenarrative"),
#'                     parser=nGram(2, ignoreCase=TRUE, 
#'                                  punctuation="[-.,?\\!:;~()]+"))
#'                     
#' # compute term-document-matrix of all 2-word combinations of Dallas police crime reports
#' # for each type of offense status
#' tdm2 = computeTfIdf(channel=NULL, tableName="public.dallaspoliceall", docId="offensestatus", 
#'                     textColumns=c("offensedescription", "offensenarrative", "offenseweather"),
#'                     parser=token(2), 
#'                     where="offensestatus NOT IN ('System.Xml.XmlElement', 'C')")
#'                     
#' # include only top 100 ranked 2-word ngrams for each 4-digit zip into resulting 
#' # term-document-matrix using rank function  
#' tdm3 = computeTfIdf(channel=NULL, tableName="public.dallaspoliceall", 
#'                     docId="substr(offensezip, 1, 4)", 
#'                     textColumns=c("offensedescription", "offensenarrative"),
#'                     parser=nGram(2), top=100)
#'                     
#' # same but get top 10% ranked terms using percent rank function
#' tdm4 = computeTfIdf(channel=NULL, tableName="public.dallaspoliceall", 
#'                     docId="substr(offensezip, 1, 4)", 
#'                     textColumns=c("offensedescription", "offensenarrative"),
#'                     parser=nGram(1), top=0.10, rankFunction="percentrank")
#' 
#' }
computeTfIdf <- function(channel, tableName, docId, textColumns, parser, 
                         top = NULL, rankFunction = 'rank',
                         idSep = '-', idNull = '(null)',
                         adjustDocumentCount = FALSE, where = NULL, 
                         stopwords = NULL, test = FALSE) {
  
  rankFunction = match.arg(rankFunction, c('rank', 'rownumber', 'row', 'denserank', 'percentrank'))
  
  if (missing(tableName)) {
    stop("Table name must be specified.")
  }
  
  if (missing(docId)) {
    stop("Doc id must be specified.")
  }
  
  if (missing(textColumns) || length(textColumns)==0) {
    stop("Text columns must be specified.")
  }
  
  isValidConnection(channel, test)

  windowFunction = getWindowFunction(rankFunction)
  
  where_clause = makeWhereClause(where)
  
  derivedDocId = makeDocumentId(docId, idSep, idNull)
  
  # validate number of documents > 1 for TF-IDF
  # and adjust 2 to 3 if requested
  if (!test) {
    countSql = paste0("SELECT COUNT(DISTINCT(", derivedDocId, ")) count ", " FROM ", tableName, where_clause)
    docCount = toaSqlQuery(channel, countSql)$count[[1]]
    if (docCount < 2)
      stop("Can't compute TF-IDF for single document. Use 'computeTf` that computes term frequency instead.")
    
    # adjust for 2 documents 
    increaseByOne = ifelse(adjustDocumentCount && docCount == 2, " + 1 ", " ")
  }else
    increaseByOne = " "
  
  textSelectSQL = parseTextSQL(parser, tableName, derivedDocId, textColumns, where)
  
  sql = paste0(
    "SELECT * FROM 
       (SELECT *, ", windowFunction, " OVER (PARTITION BY docid ORDER BY tf_idf DESC) rank FROM TF_IDF(
         ON TF(
           ON (SELECT docid, term FROM ( ", textSelectSQL, " ) t ) PARTITION BY docid
            ) AS TF PARTITION BY term
         ON ( SELECT COUNT(DISTINCT(", derivedDocId, ")) ", increaseByOne," FROM ", tableName, where_clause, " 
            ) AS doccount dimension
       )) t2", makeRankFilter(top)
    )
  
  if (test) 
    return (sql)
  else {
    rs = toaSqlQuery(channel, sql, stringsAsFactors = FALSE, na.strings="")
  }
  
  rs = removeStopWords(rs, stopwords)
  
  x = makeSimpleTripletMatrix(rs, 'tf_idf', "ti")
  return (x)
}


makeDocumentId <- function(docId, idSep, idNull) {
  
  collapse = paste0(" || '", idSep, "' || ")
  derivedId = paste0("COALESCE(CAST(", docId, " AS varchar), '", idNull, "')", collapse = collapse)
  
}

TermDocumentMatrix_classes <-
  c("toaTermDocumentMatrix", "TermDocumentMatrix", "simple_triplet_matrix")

makeSimpleTripletMatrix <- function(result_set, weight_name, weighting = "tf") {
  
  terms = result_set$term
  docs = result_set$docid
  weights = result_set[, weight_name]
  
  allTerms = as.character(sort(unique(terms)))
  allDocs = as.character(sort(unique(docs)))
  i = match(terms, allTerms)
  j = match(docs, allDocs)
  
  m = slam::simple_triplet_matrix(i = i, j = j, v = as.numeric(weights),
                            nrow = length(allTerms),
                            ncol = length(allDocs),
                            dimnames =
                              list(Terms = allTerms,
                                   Docs = allDocs))
  
  m$rs = result_set
  class(m) <- TermDocumentMatrix_classes
  attr(m, "Weighting") <- weighting
  
  return(m)
}


makeRankFilter <- function(top) {
  
  ifelse(!is.null(top) && is.numeric(top) && top > 0, 
         paste0(" WHERE rank <= ", as.character(top)), 
         "")
}

removeStopWords <- function(rs, stopwords, ignore.case=TRUE) {
  
  if (is.null(stopwords) || length(stopwords) == 0)
    return (rs)
  
  if (!ignore.case) {
    stopwords = c(tolower(stopwords), toupper(stopwords))
  }
  
  result = rs 
  for(sw in unique(stopwords)) {
    idx = grep(paste0("\\b",sw,"\\b"), result$term, ignore.case=ignore.case)
    if (length(idx)>0) result = result[-idx, ]
  }
  
  return (result)
}