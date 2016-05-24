## ------------------------------------------------------------------------
library(textreuse)
dir <- system.file("extdata/ats", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir, tokenizer = tokenize_ngrams, n = 5,
                          progress = FALSE)

## ------------------------------------------------------------------------
jaccard_similarity(corpus[["remember00palm"]], 
                   corpus[["remembermeorholy00palm"]])

## ----eval=FALSE----------------------------------------------------------
#  comparisons <- pairwise_compare(corpus, jaccard_similarity, progress = FALSE)
#  comparisons[1:4, 1:4]

## ---- echo=FALSE---------------------------------------------------------
comparisons <- pairwise_compare(corpus, jaccard_similarity, progress = FALSE)
round(comparisons[1:3, 1:3], digits = 3)

## ------------------------------------------------------------------------
candidates <- pairwise_candidates(comparisons)
candidates[candidates$score > 0.1, ]

## ----eval=FALSE----------------------------------------------------------
#  vignette("minhash", package = "textreuse")

