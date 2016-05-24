### R code from vignette source 'backtestGraphics.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: backtestGraphics.Rnw:93-103
###################################################
library(backtestGraphics)
library(dplyr)

## use the following dplyr code to show a more interesting section of
## the data frame

equity %>% filter(name %in% unique(equity$name)[1:3],
                  date %in% as.Date(c("2005-05-02", "2005-05-03",
                                      "2005-05-04"))) %>%
      arrange(date, name)


###################################################
### code chunk number 2: backtestGraphics.Rnw:123-130
###################################################

## use the following dplyr code to show a more interesting section of
## the data frame

(commodity %>% filter(id   %in% c("CO", "FC", "W"),
                      date %in% as.Date(c("2003-02-03", "2003-02-07"))) %>%
       arrange(id, name))[c(1,6,15,16,17),]


###################################################
### code chunk number 3: backtestGraphics.Rnw:151-154
###################################################
credit %>% filter(name %in% unique(credit$name)[1:3],
                  date %in% unique(credit$date)[1:3]) %>%
      arrange(date, name)


