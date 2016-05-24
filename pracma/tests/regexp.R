##
##  r e g e x p . R  Test suite
##


regexp <- pracma::regexp
regexpi <- pracma::regexpi
regexprep <- pracma::regexprep
refindall <- pracma::refindall

s <- "bat cat can car COAT court cut ct CAT-scan"
pat <-  'c[aeiou]+t'
identical(regexp(s, pat)$match,
          c("cat", "cut"))
identical(regexpi(s, pat)$match,
          c("cat", "COAT", "cut", "CAT"))
identical(regexp(s, pat, once = TRUE)$match,
          c("cat"))
identical(regexp(s, pat, ignorecase = TRUE, split = TRUE)$split,
          c("bat ", " can car ", " court ", " ct ", "-scan"))

identical(regexprep(s, pat, '---'),
          c("bat --- can car COAT court --- ct CAT-scan"))
identical(regexprep(s, pat, '---', once = TRUE),
          c("bat --- can car COAT court cut ct CAT-scan"))
identical(regexprep(s, pat, '---', ignorecase = TRUE),
          c("bat --- can car --- court --- ct ----scan"))

identical(refindall("AbababaBa", 'aba'), c(3, 5))
identical(refindall("AbababaBa", 'aba', ignorecase = TRUE), c(1, 3, 5, 7))
