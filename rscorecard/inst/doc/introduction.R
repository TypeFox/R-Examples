## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(purl = NOT_CRAN)

## ---- eval = NOT_CRAN----------------------------------------------------
library(rscorecard)
df <- sc_init() %>% 
    sc_filter(region == 2, ccbasic == c(21,22,23), locale == 41:43) %>% 
    sc_select(unitid, instnm, stabbr) %>% 
    sc_year(2013) %>% 
    sc_get()
df

## ---- eval = FALSE-------------------------------------------------------
#  # NB: You must use a real key, of course...
#  sc_key('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')

## ---- eval = FALSE-------------------------------------------------------
#  # NB: You must use a real key, of course...
#  DATAGOV_API_KEY=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

## ---- eval = NOT_CRAN----------------------------------------------------
## public schools within 50 miles of midtown Nashville, TN
df <- sc_init() %>% 
    sc_filter(control == 1) %>% 
    sc_select(unitid, instnm, stabbr) %>% 
    sc_year(2013) %>% 
    sc_zip(37203, 50) %>%
    sc_get()
df

## ---- eval = NOT_CRAN----------------------------------------------------
## median earnings for students who first enrolled in a public
## college in the New England or Mid-Atlantic regions: 10 years later
df <- sc_init() %>% 
    sc_filter(control == 1, region == 1:2, ccbasic == 1:24) %>% 
    sc_select(unitid, instnm, md_earn_wne_p10) %>% 
    sc_year(2009) %>%
    sc_get()
df

