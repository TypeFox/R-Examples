## ---- echo=FALSE, message=FALSE------------------------------------------
library("dplyr")

## ------------------------------------------------------------------------
library(textreuse)
minhash <- minhash_generator(n = 240, seed = 3552)
head(minhash(c("turn tokens into", "tokens into hashes", "into hashes fast")))

## ------------------------------------------------------------------------
dir <- system.file("extdata/ats", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir, tokenizer = tokenize_ngrams, n = 5,
                          minhash_func = minhash, keep_tokens = TRUE,
                          progress = FALSE)

## ------------------------------------------------------------------------
head(minhashes(corpus[[1]]))
length(minhashes(corpus[[1]]))

## ------------------------------------------------------------------------
lsh_threshold(h = 200, b = 50)
lsh_threshold(h = 240, b = 80)

## ------------------------------------------------------------------------
lsh_probability(h = 240, b = 80, s = 0.25)
lsh_probability(h = 240, b =  80, s = 0.75)

## ------------------------------------------------------------------------
buckets <- lsh(corpus, bands = 80, progress = FALSE)
buckets

## ------------------------------------------------------------------------
baxter_matches <- lsh_query(buckets, "calltounconv00baxt")
baxter_matches
candidates <- lsh_candidates(buckets)
candidates

## ------------------------------------------------------------------------
lsh_compare(candidates, corpus, jaccard_similarity, progress = FALSE)

