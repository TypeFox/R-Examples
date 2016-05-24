context("computeTf")

test_that("computeTf throws errors", {
  
  expect_error(computeTf(NULL, "table1", "id", c('text1','text2'), nGram(3), weighting="NO SUCH WEIGHTING"),
               "'arg' should be one of \"raw\", \"bool\", \"binary\", \"log\", \"augment\", \"normal\"")
  
  expect_error(computeTf(NULL, "table1", "id", c('text1','text2'), nGram(3), weighting="normal",
                         rankFunction = "NO SUCH RANK FUNCTION"),
               "'arg' should be one of \"rank\", \"rownumber\", \"row\", \"denserank\", \"percentrank\"")
  
  expect_error(computeTf(NULL), "Table name must be specified.")
  
  expect_error(computeTf(NULL, "table1"), "Doc id must be specified.")
  
  expect_error(computeTf(NULL, "table1", "id"), "Text columns must be specified.")
  
  expect_error(computeTf(NULL, "table1", "id", character(0)), "Text columns must be specified.")
  
  expect_error(computeTf(NULL, "table1", "id", c("text1","text2")), "Connection is not valid RODBC object.")
  
})

test_that("computeTf SQL with nGram parser is correct", {
  
  expect_equal_normalized(
    computeTf(channel=NULL, tableName="public.dallaspoliceall", docId="offensestatus", 
              textColumns=c("offensedescription", "offensenarrative"),
              parser=nGram(2), where="offensestatus NOT IN ('System.Xml.XmlElement', 'C')", test=TRUE),
    "SELECT * FROM (SELECT *, RANK() OVER (PARTITION BY docid ORDER BY tf DESC) rank FROM TF(
         ON (SELECT docid, term 
               FROM ( SELECT __doc_id__ docid, ngram as term, frequency 
                    FROM nGram (
                      ON (SELECT COALESCE(CAST(offensestatus AS varchar), '(null)') __doc_id__, 
                                 offensedescription || ' ' || offensenarrative __text_column__  
                            FROM public.dallaspoliceall 
                           WHERE offensestatus NOT IN ('System.Xml.XmlElement', 'C') )      
                      TEXT_COLUMN('__text_column__') 
                      DELIMITER('[ \\t\\b\\f\\r]+')
                      GRAMS(2) 
                      OVERLAPPING('true')
                      CASE_INSENSITIVE('false')        
                      ACCUMULATE('__doc_id__')
                    ) 
              WHERE length(ngram) >= 3
               ) t
            ) PARTITION BY docid
         FORMULA('normal')
     )) t2"
    )
  
  expect_equal_normalized(
    computeTf(channel=NULL, tableName="public.dallaspoliceall", docId="substr(offensezip, 1, 4)",
              textColumns=c("offensedescription", "offensenarrative"),
              parser=nGram(2, minLength=2), test=TRUE),
    "SELECT * FROM (SELECT *, RANK() OVER (PARTITION BY docid ORDER BY tf DESC) rank FROM TF(
         ON (SELECT docid, term
               FROM ( SELECT __doc_id__ docid, ngram as term, frequency
                    FROM nGram (
                      ON (SELECT COALESCE(CAST(substr(offensezip, 1, 4) AS varchar), '(null)') __doc_id__, offensedescription || ' ' || offensenarrative __text_column__ 
                            FROM public.dallaspoliceall )
                      TEXT_COLUMN('__text_column__')
                      DELIMITER('[ \\t\\b\\f\\r]+')
                      GRAMS(2)
                      OVERLAPPING('true')
                      CASE_INSENSITIVE('false')
                      ACCUMULATE('__doc_id__')
                    ) 
              WHERE length(ngram) >= 5
                ) t
              ) PARTITION BY docid
            FORMULA('normal')
      )) t2"
    )
})

test_that("computeTf SQL with token parser is correct", {
  
  expect_equal_normalized(
    computeTf(channel=NULL, tableName="public.dallaspoliceall", docId="(extract('hour' from offensestarttime)/6)::int%4",
              textColumns=c("offensedescription", "offensenarrative"),
              parser=token(1, punctuation="[-\\\\\\[.,?\\!:;~()\\\\\\]]+", stopWords=TRUE),
              where="offensenarrative IS NOT NULL", test=TRUE),
    "SELECT * FROM (SELECT *, RANK() OVER (PARTITION BY docid ORDER BY tf DESC) rank FROM TF(
         ON (SELECT docid, term 
               FROM ( 
                 WITH tokens AS
                   (SELECT *
                      FROM text_parser(
                        ON (SELECT COALESCE(CAST((extract('hour' from offensestarttime)/6)::int%4 AS varchar), '(null)') __doc_id__, 
                                   offensedescription || ' ' || offensenarrative __text_column__ 
                              FROM public.dallaspoliceall WHERE offensenarrative IS NOT NULL )
                        PARTITION BY __doc_id__ 
                        TEXT_COLUMN('__text_column__')
                        DELIMITER('[ \\t\\b\\f\\r]+')
                        CASE_INSENSITIVE('false')
                        STEMMING('false')
                        PUNCTUATION('[-\\\\\\[.,?\\!:;~()\\\\\\]]+')
                        REMOVE_STOP_WORDS('true')
                        ACCUMULATE('__doc_id__')
                        TOTAL('false')
                        LIST_POSITIONS('false')
                        TOKEN_COLUMN_NAME('term')
                        FREQUENCY_COLUMN_NAME('frequency')
                        OUTPUT_BY_WORD('true')
                      )
                     WHERE length(term) >= 1)
                 SELECT t1.__doc_id__ docid, term term, COUNT(*) frequency
                   FROM tokens t1
                  GROUP BY 1, 2
              ) t
            ) 
         PARTITION BY docid
         FORMULA('normal')
      )) t2"
    )
  
  expect_equal_normalized(
    computeTf(channel=NULL, tableName="public.dallaspoliceall", docId="offensebeat",
              textColumns=c("offensedescription", "offensenarrative", "offenseweather"),
              parser=token(1, stopWords="english.dat"), where="offensenarrative IS NOT NULL AND offenseweather IS NOT NULL",
              test=TRUE),
    "SELECT * FROM (SELECT *, RANK() OVER (PARTITION BY docid ORDER BY tf DESC) rank FROM TF(
         ON (SELECT docid, term 
               FROM ( 
                 WITH tokens AS
                   (SELECT *
                      FROM text_parser(
                        ON (SELECT COALESCE(CAST(offensebeat AS varchar), '(null)') __doc_id__, 
                                   offensedescription || ' ' || offensenarrative || ' ' || offenseweather __text_column__ 
                              FROM public.dallaspoliceall WHERE offensenarrative IS NOT NULL AND offenseweather IS NOT NULL )
                        PARTITION BY __doc_id__ 
                        TEXT_COLUMN('__text_column__')
                        DELIMITER('[ \\t\\b\\f\\r]+')
                        CASE_INSENSITIVE('false')
                        STEMMING('false')
                        REMOVE_STOP_WORDS('true')
                        STOP_WORDS('english.dat')
                        ACCUMULATE('__doc_id__')
                        TOTAL('false')
                        LIST_POSITIONS('false')
                        TOKEN_COLUMN_NAME('term')
                        FREQUENCY_COLUMN_NAME('frequency')
                        OUTPUT_BY_WORD('true')
                      )
                     WHERE length(term) >= 1)
                 SELECT t1.__doc_id__ docid, term term, COUNT(*) frequency
                   FROM tokens t1
                  GROUP BY 1, 2
              ) t
            ) 
         PARTITION BY docid
         FORMULA('normal')
      )) t2"
  )
  
  expect_equal_normalized(
    computeTf(channel=NULL, tableName="public.dallaspoliceall", docId="offensebeat",
              textColumns=c("offensedescription", "offensenarrative"),
              weighting='bool',
              parser=token(3, tokenSep="-+-", ignoreCase=TRUE, delimiter="[ \\t]+",
                           punctuation="[-\\\\\\[.,?\\!:;~()\\\\\\]]+", stemming=TRUE, 
                           stopWords=TRUE, minLength=4), 
              where="offensenarrative IS NOT NULL", test=TRUE),
    "SELECT * FROM (SELECT *, RANK() OVER (PARTITION BY docid ORDER BY tf DESC) rank FROM TF(
         ON (SELECT docid, term
               FROM (
                 WITH tokens AS
                   (SELECT *
                      FROM text_parser(
                        ON (SELECT COALESCE(CAST(offensebeat AS varchar), '(null)') __doc_id__, 
                                   offensedescription || ' ' || offensenarrative __text_column__
                              FROM public.dallaspoliceall WHERE offensenarrative IS NOT NULL )
                        PARTITION BY __doc_id__
                        TEXT_COLUMN('__text_column__')
                        DELIMITER('[ \\t]+')
                        CASE_INSENSITIVE('true')
                        STEMMING('true')
                        PUNCTUATION('[-\\\\\\[.,?\\!:;~()\\\\\\]]+')
                        REMOVE_STOP_WORDS('true')
                        ACCUMULATE('__doc_id__')
                        TOTAL('false')
                        LIST_POSITIONS('false')
                        TOKEN_COLUMN_NAME('term')
                        FREQUENCY_COLUMN_NAME('frequency')
                        OUTPUT_BY_WORD('true')
                      )
                     WHERE length(term) >= 4)
                   SELECT t1.__doc_id__ docid, t1.term || '-+-' || t2.term || '-+-' || t3.term term, COUNT(*) frequency
                  FROM tokens t1 JOIN 
                       tokens t2 ON (t1.__doc_id__ = t2.__doc_id__ AND t1.term < t2.term) JOIN
                       tokens t3 ON (t2.__doc_id__ = t3.__doc_id__ AND t2.term < t3.term)
                 GROUP BY 1, 2
              ) t
            )
         PARTITION BY docid
         FORMULA('bool')
     )) t2"
  )
})

test_that("computeTf SQL with top ranking is correct", {
  
  expect_equal_normalized(
    computeTf(channel=NULL, tableName="public.dallaspoliceall", docId="offensestatus", 
              textColumns=c("offensedescription", "offensenarrative"),
              parser=nGram(2), top=100, rankFunction="denserank",
              where="offensestatus NOT IN ('System.Xml.XmlElement', 'C')", test=TRUE),
    "SELECT * FROM (SELECT *, DENSE_RANK() OVER (PARTITION BY docid ORDER BY tf DESC) rank FROM TF(
         ON (SELECT docid, term 
               FROM ( SELECT __doc_id__ docid, ngram as term, frequency 
                    FROM nGram (
                      ON (SELECT COALESCE(CAST(offensestatus AS varchar), '(null)') __doc_id__, 
                                 offensedescription || ' ' || offensenarrative __text_column__  
                            FROM public.dallaspoliceall 
                           WHERE offensestatus NOT IN ('System.Xml.XmlElement', 'C') )      
                      TEXT_COLUMN('__text_column__') 
                      DELIMITER('[ \\t\\b\\f\\r]+')
                      GRAMS(2) 
                      OVERLAPPING('true')
                      CASE_INSENSITIVE('false')        
                      ACCUMULATE('__doc_id__')
                    ) 
              WHERE length(ngram) >= 3
               ) t
            ) PARTITION BY docid
         FORMULA('normal')
     )) t2
     WHERE rank <= 100"
  )
  
  expect_equal_normalized(
    computeTf(channel=NULL, tableName="public.dallaspoliceall", docId="offensestatus", 
              textColumns=c("offensedescription", "offensenarrative"),
              parser=nGram(2), top=10, 
              where="offensestatus NOT IN ('System.Xml.XmlElement', 'C')", test=TRUE),
    "SELECT * FROM (SELECT *, RANK() OVER (PARTITION BY docid ORDER BY tf DESC) rank FROM TF(
         ON (SELECT docid, term 
               FROM ( SELECT __doc_id__ docid, ngram as term, frequency 
                    FROM nGram (
                      ON (SELECT COALESCE(CAST(offensestatus AS varchar), '(null)') __doc_id__, 
                                 offensedescription || ' ' || offensenarrative __text_column__  
                            FROM public.dallaspoliceall 
                           WHERE offensestatus NOT IN ('System.Xml.XmlElement', 'C') )      
                      TEXT_COLUMN('__text_column__') 
                      DELIMITER('[ \\t\\b\\f\\r]+')
                      GRAMS(2) 
                      OVERLAPPING('true')
                      CASE_INSENSITIVE('false')        
                      ACCUMULATE('__doc_id__')
                    ) 
              WHERE length(ngram) >= 3
               ) t
            ) PARTITION BY docid
         FORMULA('normal')
     )) t2
     WHERE rank <= 10"
  )
  
  expect_equal_normalized(
    computeTf(channel=NULL, tableName="public.dallaspoliceall", docId="offensestatus", 
              textColumns=c("offensedescription", "offensenarrative"),
              parser=nGram(2), top=0.10, rankFunction="percentrank",
              where="offensestatus NOT IN ('System.Xml.XmlElement', 'C')", test=TRUE),
    "SELECT * FROM (SELECT *, PERCENT_RANK() OVER (PARTITION BY docid ORDER BY tf DESC) rank FROM TF(
         ON (SELECT docid, term 
               FROM ( SELECT __doc_id__ docid, ngram as term, frequency 
                    FROM nGram (
                      ON (SELECT COALESCE(CAST(offensestatus AS varchar), '(null)') __doc_id__, 
                                 offensedescription || ' ' || offensenarrative __text_column__  
                            FROM public.dallaspoliceall 
                           WHERE offensestatus NOT IN ('System.Xml.XmlElement', 'C') )      
                      TEXT_COLUMN('__text_column__') 
                      DELIMITER('[ \\t\\b\\f\\r]+')
                      GRAMS(2) 
                      OVERLAPPING('true')
                      CASE_INSENSITIVE('false')        
                      ACCUMULATE('__doc_id__')
                    ) 
              WHERE length(ngram) >= 3
               ) t
            ) PARTITION BY docid
         FORMULA('normal')
     )) t2
     WHERE rank <= 0.1"
  )
  
})

test_that("stop words in R works", {
  rs = data.frame(term=c("a","b","c", "the", "THE"), freq=c(1,2,3,4,5), stringsAsFactors=FALSE)
  expect_equal(removeStopWords(rs, NULL), rs)
  
  expect_equal(removeStopWords(rs, character(0)), rs)
  
  expect_equal(removeStopWords(rs, c("a","b","c", "the")), data.frame(term=character(0), freq=numeric(0), stringsAsFactors=FALSE))
  
  expect_equal(removeStopWords(rs, c("a","b", "the")), rs[3,])
  
  expect_equal(removeStopWords(rs, c("a","b", "the"), ignore.case=FALSE), rs[3,])
})