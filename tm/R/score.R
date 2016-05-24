tm_term_score <-
function(x, terms, FUN)
    UseMethod("tm_term_score", x)

tm_term_score.term_frequency <-
function(x, terms, FUN = function(x) sum(x, na.rm = TRUE))
    FUN(x[match(terms, names(x), nomatch = 0L)])

tm_term_score.PlainTextDocument <-
function(x, terms, FUN = function(x) sum(x, na.rm = TRUE))
    tm_term_score(termFreq(x, control = list(tolower = FALSE,
                                             removePunctuation = TRUE,
                                             wordLengths = c(1, Inf))),
                  terms, FUN)

tm_term_score.TermDocumentMatrix <-
function(x, terms, FUN = slam::col_sums)
    FUN(x[match(terms, Terms(x), nomatch = 0L), ])

tm_term_score.DocumentTermMatrix <-
function(x, terms, FUN = slam::row_sums)
    FUN(x[, match(terms, Terms(x), nomatch = 0L)])
