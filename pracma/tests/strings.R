##
##  s t r i n g s . R  tests
##

strcat <- pracma::strcat
strcmp <- pracma::strcmp
strcmpi <- pracma::strcmpi

strcmp(" empty", " empty")               # T
!strcmp("empty ", "empty")               # F
!strcmp("foobar", "barfoo")              # F
!strcmp("string", "String")              # F
!strcmp(c("yes", "no"), c("yes", "on"))  # F
!strcmp(c("abc", "abc"), c("abc"))       # F
strcmp(c("yes", "no"), c("yes", "no"))   # T

strcmpi("string", "String")              # T
strcmpi(c("yes", "no"), c("Yes", "No"))  # T

blanks <- pracma::blanks
deblank <- pracma::deblank
strTrim <- pracma::strTrim
strjust <- pracma::strjust
strRep <- pracma::strRep

identical(c(blanks(0), blanks(1), blanks(2)), c("", " ", "  "))
s <- c("  abc", "abc   ", " abc ", " a b c ", "abc", "a b c")
identical(deblank(s), c("  abc", "abc", " abc", " a b c", "abc", "a b c"))
identical(strTrim(s), c("abc", "abc", "abc", "a b c", "abc", "a b c"))
identical(strjust(s, justify = "center"),
          c(" abc ", " abc ", " abc ", "a b c", " abc ", "a b c"))
s <- c('This is a good example.', "He has a good character.",
       'This is good, good food.', "How goodgood this is!")
identical(strRep(s, 'good', 'great'),
          c('This is a great example.', "He has a great character.",
            'This is great, great food.', "How greatgreat this is!"))
