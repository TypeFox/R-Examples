library(microbenchmark)
library(text2vec)
data("movie_review")
tokens <- movie_review$review %>% tolower %>% word_tokenizer
times = 5
microbenchmark(v1 <- create_vocabulary(itoken(tokens), ngram = c(1, 1)),  times = times)
# 0.406   0.008   0.415
# 0.391   0.005   0.392
microbenchmark(v2 <- create_vocabulary(itoken(tokens), ngram = c(1, 2)), times = times )
# 1.845   0.030   1.872
# 1.416   0.018   1.430

# system.time(dtm <- create_dtm(itoken(tokens), vocab_vectorizer(v1)) )
microbenchmark(corpus <- create_corpus(itoken(tokens), vocab_vectorizer(v1)), times = times )
# 0.504   0.008   0.508
microbenchmark(corpus <- create_corpus(itoken(tokens), vocab_vectorizer(v2)), times = times )
# 1.885   0.044   1.926


microbenchmark(corpus <- create_corpus(itoken(tokens), hash_vectorizer(ngram = c(1, 1))),  times = times )
# 0.419   0.007   0.423

microbenchmark(corpus <- create_corpus(itoken(tokens), hash_vectorizer(ngram = c(1, 2))),  times = times )
# 1.290   0.031   1.317
dtm <- corpus$get_dtm()
