## ----message=FALSE-------------------------------------------------------
library(htmltab)

url <- "http://en.wikipedia.org/wiki/Demography_of_the_United_Kingdom"
ukLang <- htmltab(doc = url, which = "//th[text() = 'Ability']/ancestor::table")
head(ukLang)

## ----message=FALSE-------------------------------------------------------
library(tidyr)
library(magrittr)

ukLang %<>% gather(key, value, -Ability)

## ------------------------------------------------------------------------
ukLang %>% separate(key, into = c("region", "language", "statistic"), sep = " >> ") %>% head

## ------------------------------------------------------------------------
url <- "http://en.wikipedia.org/wiki/New_Zealand_general_election,_2002"
xp <- "//caption[starts-with(text(), 'Electorate results')]/ancestor::table"

nz1 <- htmltab(doc = url, which = xp, encoding = "UTF-8")
head(nz1)

## ------------------------------------------------------------------------
nz2 <- htmltab(doc = url, which = xp, header = 1 + "//tr/td[@colspan = '10']",
               body = "//tr[./td[not(@colspan = '10')]]", encoding = "UTF-8")

## ------------------------------------------------------------------------
tail(nz2, 9)

