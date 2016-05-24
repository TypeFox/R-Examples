## user-specifiable lookup-functions
match_ignorecase <-
function(lookup, entry, ...)
    tolower(lookup) %in% tolower(entry)

match_exact <-
function(lookup, entry, ...)
    lookup %in% entry

match_partial <-
function(lookup, entry, ...)
    !is.na(pmatch(lookup, entry, ...))

match_partial_ignorecase <-
function(lookup, entry, ...)
    !is.na(pmatch(tolower(lookup), tolower(entry), ...))

match_regexp <-
function(lookup, entry, ...)
    length(grep(lookup, entry, ...)) > 0
