## ----init, echo=FALSE, message=FALSE-------------------------------------
library(dplyr, quietly = T)
library(cbsodataR)

## ----get_table_list, message=FALSE---------------------------------------
tables <- get_table_list(Language="en") # retrieve only enlgish tables

tables %>% 
  select(Identifier, ShortTitle) %>% 
  head 

## ----get_meta, message=FALSE---------------------------------------------
m <- get_meta('71509ENG')
m

## ---- meta2--------------------------------------------------------------
names(m)

## ---- get_data2, message=FALSE-------------------------------------------
get_data('71509ENG') %>% 
  select(2:5) %>%  # select column 2 to 5 (for demonstration purpose)
  head

## ---- get_data, message=FALSE--------------------------------------------
get_data('71509ENG', recode = FALSE) %>% 
  select(2:5) %>% 
  head

## ---- get_data3, message=FALSE-------------------------------------------
  get_data('71509ENG', Periods='2000JJ00') %>% 
  select(2:5) %>% 
  head

