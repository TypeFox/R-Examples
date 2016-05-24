# Author: Ingo Feinerer
# Filters

tm_filter <-
function(x, FUN, ...)
    UseMethod("tm_filter", x)
tm_filter.PCorpus <- tm_filter.VCorpus <-
function(x, FUN, ...)
    x[tm_index(x, FUN, ...)]

tm_index <-
function(x, FUN, ...)
    UseMethod("tm_index", x)
tm_index.PCorpus <- tm_index.VCorpus <-
function(x, FUN, ...)
    unlist(mclapply(content(x), FUN, ...))
