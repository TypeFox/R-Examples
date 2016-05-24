## ----setup, echo = FALSE-------------------------------------------------
library(knitr)
opts_chunk$set(cache = TRUE, message = FALSE)

## ------------------------------------------------------------------------
library(dplyr)
library(fuzzyjoin)
data(misspellings)

misspellings

## ----words---------------------------------------------------------------
# use the dictionary of words from the qdapDictionaries package,
# which is based on the Nettalk corpus.
library(qdapDictionaries)
words <- tbl_df(DICTIONARY)

words

## ----sub_misspellings----------------------------------------------------
set.seed(2016)
sub_misspellings <- misspellings %>%
  sample_n(1000)

## ----joined, dependson = c("words", "sub_misspellings")------------------
joined <- sub_misspellings %>%
  stringdist_inner_join(words, by = c(misspelling = "word"), max_dist = 1)

## ----dependson = "joined"------------------------------------------------
joined

## ----dependson = "joined"------------------------------------------------
joined %>%
  count(misspelling, correct)

## ----dependson = "joined"------------------------------------------------
which_correct <- joined %>%
  group_by(misspelling, correct) %>%
  summarize(guesses = n(), one_correct = any(correct == word))

which_correct

# percentage of guesses getting at least one right
mean(which_correct$one_correct)

# number uniquely correct (out of the original 1000)
sum(which_correct$guesses == 1 & which_correct$one_correct)

## ----left_joined, dependson = "misspellings"-----------------------------
left_joined <- sub_misspellings %>%
  stringdist_left_join(words, by = c(misspelling = "word"), max_dist = 1)

left_joined

left_joined %>%
  filter(is.na(word))

## ----left_joined2, dependson = "misspellings"----------------------------
left_joined2 <- sub_misspellings %>%
  stringdist_left_join(words, by = c(misspelling = "word"), max_dist = 2)

left_joined2

left_joined2 %>%
  filter(is.na(word))

