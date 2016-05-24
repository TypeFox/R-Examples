getKeywords <- function (pos = 2:length(search()))
{
	## Get a sorted list of unique function names for libraries loaded
	## in positions provided by pos
	res <- NULL
	for (i in pos) {
		if (search()[i] == "package:base")  # Use builtins() instead
			res <- c(res, builtins()) else
			res <- c(res, as.character(getFunctions(i)))
	}
	## Sort res and return only unique names
	res <- sort(res[!duplicated(res)])
	## Eliminate items containing <-, __, -, !, $, %, &, |, *, +, /, :, [ or =
	searchit <- c("<-", "__", "-", "!", "[$]", "[%]", "[&]", "[|]", "[*]",
		"[+]", "[/]", ":", "[[]", "=")
	for (i in 1:length(searchit)) {
		elim <- grep(searchit[i], res)
		if (length(elim) > 0) res <- res[-elim]
	}
	## Eliminate some other items (reserved keywords already introduced in
	## keyword1 list, and other stuff)
	reserved <- c("break", "else", "FALSE", "for", "function", "if", "in",
		"Inf", "NA", "NaN", "next", "NULL", "repeat", "TRUE", "while", "(",
		"?", "@", "^", "{", "~", "<", ">")
	for (i in 1:length(reserved))
		res <- res[res != reserved[i]]
	return(res)
}
