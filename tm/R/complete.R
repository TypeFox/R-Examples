# Author: Ingo Feinerer

stemCompletion <-
function(x, dictionary,
         type = c("prevalent", "first", "longest",
                  "none", "random", "shortest"))
{
    if (inherits(dictionary, "Corpus"))
        dictionary <- unique(unlist(lapply(dictionary, words)))

    type <- match.arg(type)
    possibleCompletions <- lapply(x, function(w) grep(sprintf("^%s", w),
                                                      dictionary,
                                                      value = TRUE))
    switch(type,
           first = {
               setNames(sapply(possibleCompletions, "[", 1), x)
           },
           longest = {
               ordering <-
                   lapply(possibleCompletions,
                          function(x) order(nchar(x), decreasing = TRUE))
               possibleCompletions <-
                   mapply(function(x, id) x[id], possibleCompletions,
                          ordering, SIMPLIFY = FALSE)
               setNames(sapply(possibleCompletions, "[", 1), x)
           },
           none = {
               setNames(x, x)
           },
           prevalent = {
               possibleCompletions <-
                   lapply(possibleCompletions,
                          function(x) sort(table(x), decreasing = TRUE))
               n <- names(sapply(possibleCompletions, "[", 1))
               setNames(if (length(n)) n else rep(NA, length(x)), x)
           },
           random = {
               setNames(sapply(possibleCompletions, function(x) {
                   if (length(x)) sample(x, 1) else NA
               }), x)
           },
           shortest = {
               ordering <- lapply(possibleCompletions, function(x) order(nchar(x)))
               possibleCompletions <-
                   mapply(function(x, id) x[id], possibleCompletions,
                          ordering, SIMPLIFY = FALSE)
               setNames(sapply(possibleCompletions, "[", 1), x)
           }
           )
}
