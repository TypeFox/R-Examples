context('str_collapse_whitespace')

"A  B"    %>% str_collapse_whitespace %>% expect_equal("A B")
"A__B__C" %>% str_collapse_whitespace %>% expect_equal("A_B_C")
"A  B__C" %>% str_collapse_whitespace %>% expect_equal("A B_C")
"A _B_ C" %>% str_collapse_whitespace %>% expect_equal("A _B_ C")

"A _B_ C" %>% str_collapse_whitespace('[\\s-_]') %>% expect_equal("A_B C")
"A _B_ C" %>% str_collapse_whitespace() %>% expect_equal("A _B_ C")


"A _B_ C" %>% 
  str_collapse_whitespace(c("\\s", "_")) %>%   # No match
  expect_equal("A _B_ C")   

context('..using regex class for pattern')
"A _B_ C" %>% 
  str_collapse_whitespace('[\\s-_]') %>%   # No match
  expect_equal("A_B C")  
