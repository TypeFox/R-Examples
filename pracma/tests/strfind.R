##
##  s t r f i n d . R  Test suite
##


strfind <- pracma::strfind
strfindi <- pracma::strfindi
#findstr <- pracma::findstr

identical(strfind("", "aba"), NULL)
identical(strfind("ab", "aba"), NULL)
identical(strfind("aba", "aba"), 1)
identical(strfind("ababa", "aba"), c(1, 3))
identical(strfind("ababa", "aba", overlap=FALSE), 1)

identical(strfindi("ABA", "aba"), 1)
identical(strfindi("aba", "ABA"), 1)
identical(strfindi("ABABA", "aba"), c(1, 3))
identical(strfindi("aBaBa", "AbA", overlap=FALSE), 1)
