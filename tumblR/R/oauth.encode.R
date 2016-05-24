oauth.encode <-
function(x) vapply(x, oauth.encode1, character(1))
