WN_filter_types <-
    c("ContainsFilter", "EndsWithFilter", "ExactMatchFilter",
      "RegexFilter", "SoundFilter", "StartsWithFilter",
      "WildcardFilter")

getFilterTypes <-
function()
    WN_filter_types

getTermFilter <-
function(type, word, ignoreCase)
{
    type <- .expand_filter_type(type[1L])
    .jnew(paste("com.nexagis.jawbone.filter", type, sep = "."),
          word, ignoreCase)
}

.expand_filter_type <-
function(x)
{
    y <- pmatch(tolower(x), tolower(WN_filter_types))
    if(is.na(y))
        stop(sprintf("Unknown filter type '%s'", x))
    WN_filter_types[y]
}
