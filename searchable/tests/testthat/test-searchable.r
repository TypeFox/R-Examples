library(testthat)
library(searchable)
library(magrittr)


# CREATE UTF map of A-z
  utf8 <- 65:122
  names(utf8) <- 
    intToUtf8(utf8) %>% strsplit('')  %>% extract2(1)  
  
  utf8. <- searchable(utf8)  

context( 'utf8-Extract' )
  utf8.[ "a" ]         %>% expect_equivalent(97) 
  utf8.[ fixed("a") ]  %>% expect_equivalent(97) 
  utf8.[ ignore.case("a") ] %>% expect_equivalent( c(65,97) )
  # utf8.[ reverse.lookup("65") ]

context( 'Replace' )
  utf8.[ "a" ]  <- -97      
  utf8.[ "a" ] %>% expect_equivalent(-97) 
  