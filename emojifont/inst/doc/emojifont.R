## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE---------------------------
library("ggplot2")
library("emojifont")

## ------------------------------------------------------------------------
library(emojifont)

search_emoji('smile')
emoji(search_emoji('smile'))

## ------------------------------------------------------------------------
## list available emoji fonts
list.emojifonts()

## load selected emoji font
load.emojifont('OpenSansEmoji.ttf')

## ----echo=FALSE----------------------------------------------------------
sessionInfo()

