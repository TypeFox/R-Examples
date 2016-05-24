## ----frankenstein--------------------------------------------------------
library(dplyr)
library(gutenbergr)

frankenstein <- gutenberg_download(84)

frankenstein

## ----frankenstein_dracula------------------------------------------------
frankenstein_dracula <- gutenberg_download(c(84, 345), meta_fields = "title")

frankenstein_dracula

## ------------------------------------------------------------------------
frankenstein_dracula %>%
  count(title)

## ------------------------------------------------------------------------
gutenberg_metadata

## ------------------------------------------------------------------------
gutenberg_metadata %>%
  filter(title == "Wuthering Heights")

## ------------------------------------------------------------------------
gutenberg_works()

## ------------------------------------------------------------------------
gutenberg_works(author == "Austen, Jane")

## ------------------------------------------------------------------------
gutenberg_subjects

## ------------------------------------------------------------------------
gutenberg_subjects %>%
  filter(subject == "Horror tales")

gutenberg_subjects %>%
  filter(grepl("Holmes, Sherlock", subject))

## ------------------------------------------------------------------------
gutenberg_authors

## ------------------------------------------------------------------------
library(tidytext)

words <- frankenstein_dracula %>%
  unnest_tokens(word, text)

words

word_counts <- words %>%
  anti_join(stop_words, by = "word") %>%
  count(title, word, sort = TRUE)

word_counts

