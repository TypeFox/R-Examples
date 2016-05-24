as.farsi <- function (x, ...)
UseMethod("farsi")

farsi <- function (x, ...)
UseMethod("farsi")

farsi.default <- function (x, ...)
    stop("no method defined for objects of class ",
         paste(dQuote(class(x)), collapse = ", "))

farsi.farsi <- function (x, ...)
structure(x, class = "farsi")

farsi.numeric <- function (x, ...)
structure(x, class = "farsi")

Ops.farsi<- function (e1, e2) {
  e1 <- unclass(e1)
  if (!missing(e2))
      e2 <- unclass(e2)
  structure(NextMethod(.Generic), class = "farsi")
}

as.numeric.farsi <- function(x, ...)
    unclass(x)

print.farsi <- function (x, ...) {
  print(noquote(as.character(x)))
  invisible(x)
}

sort.farsi <- function (x, decreasing = FALSE, ...) {
  structure(NextMethod("sort"), class = "farsi")
}

as.character.farsi <- local({
  helper <- local({
    check <- function (dvec)
        (d <- makeNumber(dvec)) > 0 && d < 1000
    makeNumber <- function (n)
        as.numeric(paste(n, collapse = ""))
    yekan <-
        structure(c("", "\uFBFE\uFB90", "\u062F\u0648", "\uFEB3\uFEEA", "\uFB7C\uFEEC\uFE8E\u0631", "\uFB58\uFEE7\uFE9D", "\uFEB7\uFEB5",
                    "\uFEEB\uFED3\uFE95", "\uFEEB\uFEB7\uFE95", "\uFEE7\uFEEA"),
                  .Names = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
    suffixes <- c("\uFEEB\u0632\uFE8E\u0631", "\uFEE3\uFBFE\uFEDF\uFBFE\u0648\uFEE5", "\uFEE3\uFBFE\uFEDF\uFBFE\uFE8E\u0631\u062F", "\uFEEB\u0632\uFE8E\u0631\uFEE3\uFBFE\uFEDF\uFBFE\uFE8E\u0631\u062F", "\uFEE3\uFBFE\uFEDF\uFBFE\u0648\uFEE5 \uFEE3\uFBFE\uFEDF\uFBFE\uFE8E\u0631\u062F ",
                  "\uFEE3\uFBFE\uFEDF\uFBFE\uFE8E\u0631\u062F\uFEE3\uFBFE\uFEDF\uFBFE\uFE8E\u0631\u062F", "sextillion", "septillion")
    teens <- structure(c("\u062F\uFEE9", "\uFBFE\uFE8E\u0632\u062F\uFEE9", "\u062F\u0648\uFE8E\u0632\u062F\uFEE9", "\uFEB3\uFBFE\u0632\u062F\uFEE9", "\uFB7C\uFEEC\uFE8E\u0631\u062F\uFEE9",
                         "\uFB58\uFE8E\uFEE7\u0632\u062F\uFEE9", "\uFEB7\uFE8E\uFEE7\u0632\u062F\uFEE9", "\uFEEB\uFED3\u062F\uFEE9", "\uFEEB\uFE9F\u062F\uFEE9",
                         "\uFEE7\u0648\u0632\u062F\uFEE9"),
                       .Names = c("0", "1", "2", "3", "4", "5", "6",
                       "7", "8", "9"))
    dahgan <- structure(c("\uFE91\uFBFE\uFEB3\uFE95", "\uFEB3\uFBFD ", "\uFB7C\uFEEC\uFEDD", "\uFB58\uFEE7\uFE9F\uFE8E\uFEE9", "\uFEB7\uFEBB\uFE95",
                        "\uFEEB\uFED3\uFE97\uFE8E\u062F", "\uFEEB\uFEB7\uFE97\uFE8E\u062F", "\uFEE7\u0648\u062F"),
                      .Names = c("2", "3", "4", "5", "6", "7", "8", "9"))
   sadgan <- structure(c("\uFEBB\u062F", "\u062F\u0648\uFBFE\uFEB3\uFE95", "\uFEB3\uFBFE\uFEBB\u062F", "\uFB7C\uFEEC\uFE8E\u0631\uFEBB\u062F", "\uFB58\uFE8E\uFEE7\uFEBB\u062F",
                        "\uFEB7\uFEB7\uFEBB\u062F", "\uFEEB\uFED3\uFE97\uFEBB\u062F", "\uFEEB\uFEB7\uFE97\uFEBB\u062F", "\uFEE7\uFEEC\uFEBB\u062F"),
                      .Names = c( "1", "2", "3", "4", "5", "6",
                       "7", "8", "9"))
    trim <- function (text)
        sub("^ *", "", sub(" *$", "", gsub("  +", " ", text)))

    function (x) {
      digits <- rev(strsplit(as.character(x), "")[[1]])
      nDigits <- length(digits)
      if (nDigits == 1)
          as.vector(yekan[digits])
      else if (nDigits == 2)
          if (x <= 19)
              as.vector(teens[digits[1]])
          else trim(paste(dahgan[digits[2]], 
                        if (check(digits[1])) "\u0648"
   else "",helper(as.numeric(digits[1]))))
      else if (nDigits == 3)
           trim(paste(sadgan[digits[3]],
                     if (check(digits[2]))"\u0648"
                      else "", helper(makeNumber(digits[2:1]))))
      else {
        nSuffix <- ((nDigits + 2)%/%3) - 1
        if (nSuffix > length(suffixes)) {
          warning(paste(x, "is too large!"))
          return(as.character(as.vector(x)))
        }
        trim(paste(helper(makeNumber(digits[nDigits:(3 * nSuffix + 1)])),
                   suffixes[nSuffix], if (check(digits[(3 * nSuffix):1]))
                   "\u0648" else "\u0648", helper(makeNumber(digits[(3 * nSuffix):1]))))
      }
    }

  })

  function (x, ...) {
    opts <- options(scipen = 1000)
    on.exit(options(opts))
    r <- character(length(x))
    bad <- is.na(x) | is.nan(x) | is.infinite(x)
    if (any(!bad & x%%1 != 0)) {
      warning("non-integer values rounded for display")
      x <- round(x)
    }
    if (any(n <- !bad & x < 0))
        r[n] <- paste("minus", sapply(-x[n], helper))
    if (any(z <- !bad & x == 0))
        r[z] <- "zero"
    if (any(p <- !bad & x > 0))
        r[p] <- sapply(x[p], helper)
    r[is.na(x)] <- ""
    r[is.nan(x)] <- "not a number"
    if (any(k <- x < 0 & is.infinite(x)))
        r[k] <- "minus infinity"
    if (any(k <- x > 0 & is.infinite(x)))
        r[k] <- "infinity"
    names(r) <- names(x)
    r
  }
})
