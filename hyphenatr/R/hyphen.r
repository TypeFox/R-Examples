#' Switch hyphen rules language
#'
#' See \code{References} for the source of the hyphenation rules. For valid
#' values, use \code{list_dicts()}.
#'
#' @param lang hyphen rules language
#' @references \url{https://github.com/LibreOffice/dictionaries}
#' @export
#' @note The \code{en_US} hyphenation rules dictionary is loaded by default
#' @examples
#' switch_dict("hu_HU")
switch_dict <- function(lang) {
  lang <- match.arg(lang, list_dicts())
  .pkgenv$curr_lang <- lang
  if (interactive()) message(sprintf("Switching hyphenation rules to '%s'", lang))
  invisible(.Call('hyphenatr_init', PACKAGE='hyphenatr',
                  system.file(sprintf("extdata/dicts/hyph_%s.dic", lang),
                              package="hyphenatr")))
}

#' Identify current hyphen rules language
#'
#' @export
#' @examples
#' curr_dict()
curr_dict <- function() { .pkgenv$curr_lang }

#' List available hyphenation languages rules
#'
#' The vast majority \href{https://github.com/LibreOffice/dictionaries}{OpenOffice Hunspell hyphenation rules dictionaries}
#' are named in the format \code{location.code_COUNTRY.CODE} (e.g. \code{en_US}).
#' This function returns the list of available dictionaries. Used the desired, supplied
#' result in your call to \code{\link{switch_dict}}.
#'
#' @export
#' @note The \code{en_US} hyphenation rules dictionary is loaded by default
#' @examples
#' list_dicts()
list_dicts <- function() {

  sort(gsub("(^hyph_|\\.dic$)", "",
            list.files(system.file("extdata/dicts",
                                   package="hyphenatr"))))

}

#' Hyphenate a character vector of words
#'
#' Given a character vector (one word per element), this function
#' will hyphenate the strings or return a list of separated
#' hyphenated string components.
#'
#' @param words character vector of words
#' @param simplify if \code{TRUE}, will return \code{words} with
#'        \code{=} as the hyphen character. If \code{FALSE}, will
#'        return a list of separated, hyphenated word components.
#'        If a string (e.g. "\code{-}", "\code{&shy;}" or '\code{&#173;}")
#'        will use that character for the hyphen character.)
#' @return a character vector or a list depending on the value of
#'         \code{simplify}.
#' @export
#' @note The \code{en_US} hyphenation rules dictionary is loaded by default
#' @examples
#' dat <- readLines(system.file("extdata/top10000en.txt", package="hyphenatr"))
#'
#' out1 <- hyphenate(dat)
#' out2 <- hyphenate(dat, simplify=FALSE)
#' out3 <- hyphenate(dat, simplify="-")
#' out4 <- hyphenate(dat, simplify="&shy;")
hyphenate <- function(words, simplify=TRUE) {

  out <- .Call('hyphenatr_hyphenate', words, PACKAGE = 'hyphenatr')

  if (inherits(simplify, "logical")) {
    if (!simplify) {
      stri_split_regex(out, "=")
    } else {
      out
    }
  } else {
    stri_replace_all_regex(out, "=", simplify)
  }

}

