## ----setup, echo=FALSE---------------------------------------------------
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE)

## ------------------------------------------------------------------------
library(janeaustenr)
library(dplyr)
library(stringr)

original_books <- austen_books() %>%
  group_by(book) %>%
  mutate(linenumber = row_number(),
         chapter = cumsum(str_detect(text, regex("^chapter [\\divxlc]",
                                                 ignore_case = TRUE)))) %>%
  ungroup()

original_books

## ------------------------------------------------------------------------
library(tidytext)
tidy_books <- original_books %>%
  unnest_tokens(word, text)

tidy_books

## ------------------------------------------------------------------------
data("stop_words")
tidy_books <- tidy_books %>%
  anti_join(stop_words)

## ------------------------------------------------------------------------
tidy_books %>%
  count(word, sort = TRUE) 

## ------------------------------------------------------------------------
nrcjoy <- sentiments %>%
  filter(lexicon == "nrc", sentiment == "joy")

tidy_books %>%
  filter(book == "Emma") %>%
  semi_join(nrcjoy) %>%
  count(word, sort = TRUE)

## ------------------------------------------------------------------------
library(tidyr)
bing <- sentiments %>%
  filter(lexicon == "bing") %>%
  select(-score)

janeaustensentiment <- tidy_books %>%
  inner_join(bing) %>%
  count(book, index = linenumber %/% 80, sentiment) %>%
  spread(sentiment, n, fill = 0) %>%
  mutate(sentiment = positive - negative)

## ---- fig.width=7, fig.height=7, warning=FALSE---------------------------
library(ggplot2)

ggplot(janeaustensentiment, aes(index, sentiment, fill = book)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~book, ncol = 2, scales = "free_x")

## ------------------------------------------------------------------------
bing_word_counts <- tidy_books %>%
  inner_join(bing) %>%
  count(word, sentiment, sort = TRUE) %>%
  ungroup()

bing_word_counts

## ---- fig.width=7, fig.height=6------------------------------------------
bing_word_counts %>%
  filter(n > 150) %>%
  mutate(n = ifelse(sentiment == "negative", -n, n)) %>%
  mutate(word = reorder(word, n)) %>%
  ggplot(aes(word, n, fill = sentiment)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Contribution to sentiment")

## ---- fig.height=6, fig.width=6------------------------------------------
library(wordcloud)

tidy_books %>%
  count(word) %>%
  with(wordcloud(word, n, max.words = 100))

## ----wordcloud, fig.height=5, fig.width=5--------------------------------
library(reshape2)

tidy_books %>%
  inner_join(bing) %>%
  count(word, sentiment, sort = TRUE) %>%
  acast(word ~ sentiment, value.var = "n", fill = 0) %>%
  comparison.cloud(colors = c("#F8766D", "#00BFC4"),
                   max.words = 100)

## ------------------------------------------------------------------------
austen_sentences <- austen_books() %>% 
  group_by(book) %>% 
  unnest_tokens(sentence, text, token = "sentences") %>% 
  ungroup()

## ------------------------------------------------------------------------
austen_sentences$sentence[39]

## ------------------------------------------------------------------------
bingnegative <- sentiments %>%
  filter(lexicon == "bing", sentiment == "negative")

wordcounts <- tidy_books %>%
  group_by(book, chapter) %>%
  summarize(words = n())

tidy_books %>%
  semi_join(bingnegative) %>%
  group_by(book, chapter) %>%
  summarize(negativewords = n()) %>%
  left_join(wordcounts, by = c("book", "chapter")) %>%
  mutate(ratio = negativewords/words) %>%
  filter(chapter != 0) %>%
  top_n(1)

