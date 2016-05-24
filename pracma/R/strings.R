##
##  s t r i n g s . R
##


strcat <- function(s1, s2 = NULL, collapse = "") {
	stopifnot(is.character(collapse))
	if (!is.vector(s1, mode = "character"))
	    stop("Argument 's1' must be a character vector.")

	if (is.null(s2)) {
	    paste(s1, collapse=collapse)
	} else {
	    if (!is.vector(s2, mode = "character"))
	        stop("Argument 's2' must be a character vector.")
	    else
	        paste(rep(s1, each = length(s2)), s2, sep = collapse)
    }
}

strcmp <- function(s1, s2) {
	if (!is.vector(s1, mode="character") || !is.vector(s1, mode="character"))
	    stop("Arguments 's1' and 's2' must be character vectors.")

    if (length(s1) == length(s2))
        all(s1 == s2)
    else
	    FALSE
}

strcmpi <- function(s1, s2) {
	if (!is.vector(s1, mode="character") || !is.vector(s1, mode="character"))
	    stop("Arguments 's1' and 's2' must be character vectors.")

    strcmp(tolower(s1), tolower(s2))
}

strTrim <- function(s) {
    if (! is.character(s))
        stop("Argument 's' must be a character vector.")

    sub("\\s+$", "", sub("^\\s+", "", s))
}

deblank <- function(s) {
    if (! is.character(s))
        stop("Argument 's' must be a character vector.")

    sub("\\s+$", "", s)
}

blanks <- function(n = 1) {
    stopifnot(is.numeric(n), length(n) == 1, n >= 0)
    n <- floor(n)

    paste(rep(" ", n), collapse="")
}

strjust <- function(s, justify = c("left", "right", "center")) {
    if (! is.character(s))
        stop("Argument 's' must be a character vector.")

    justify <- match.arg(justify)

    s <- strTrim(s)
    n <- length(s)
    M <- nchar(s)
    m <- max(M)

    S <- character(n)
    for (i in 1:n) {
        k <- m - M[i]
        if (justify == "left") {
            S[i] <- paste(s[i], blanks(k), sep = "", collapse="")
        } else if (justify == "right") {
            S[i] <- paste(blanks(k), s[i], sep = "", collapse="")
        } else {  # justify == "center"
            kl <- k %/% 2
            kr <- k - kl
            S[i] <- paste(blanks(kl), s[i], blanks(kr), sep = "", collapse="")
        }
   }
   return(S)
}

strRep <- function(s, old, new) {
    # Find and replace substring
    if (! is.character(s))
        stop("Argument 's' must be a character vector.")
    if (!is.character(old) || !is.character(new) ||
        length(old) != 1   || length(new) != 1)
        stop("Arguments 'old' and 'new' must be simple character strings.")

    gsub(old, new, s, fixed = TRUE)
}
