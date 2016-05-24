## ---- message = FALSE----------------------------------------------------
library(easyformatr)
library(magrittr)
library(dplyr)
library(knitr)

## ------------------------------------------------------------------------
filter_info(type == "time" & component == "base")

## ------------------------------------------------------------------------
filter_info(type == "number" & component == "base")

## ------------------------------------------------------------------------
easy_format(year, month, day, integer, octal, double)

## ------------------------------------------------------------------------
filter_info(type == "time" & component == "mutant" & base == "second")

## ------------------------------------------------------------------------
filter_info(type == "number" & component == "flag")

## ------------------------------------------------------------------------
filter_info(type == "number" & component == "option")

## ------------------------------------------------------------------------
easy_format(second %>% decimal)

## ------------------------------------------------------------------------
easy_format(list(integer, 
                 double) %>% 
              always_decimal)

## ------------------------------------------------------------------------
easy_format(list(month, 
                 list(day,
                      minute) ) %>% 
            roman)

## ------------------------------------------------------------------------
easy_format(second %>% roman) ==
  easy_format(second)

## ------------------------------------------------------------------------
easy_format(second %>% decimal) == 
  easy_format(second %>% decimal(TRUE) )

## ------------------------------------------------------------------------
easy_format(second %>% 
              decimal %>% 
              decimal(NA) ) ==
  easy_format(second)

## ------------------------------------------------------------------------
easy_format(double %>% before_decimal) == 
  easy_format(double)

## ------------------------------------------------------------------------
easy_format(double %>% before_decimal(3) %>% before_decimal) == 
  easy_format(double)

## ---- error = TRUE-------------------------------------------------------
easy_format(second %>%
              decimal %>%
              digits(1) ) 

## ------------------------------------------------------------------------
easy_format("We are the 99%")

## ------------------------------------------------------------------------
easy_format(roman("I am Spartacus"))

## ------------------------------------------------------------------------
easy_format("We", "are", "the 99%", sep = " ")

