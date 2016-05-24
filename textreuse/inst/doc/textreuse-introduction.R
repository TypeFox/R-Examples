## ----eval=FALSE----------------------------------------------------------
#  vignette("textreuse-pairwise", package = "textreuse")
#  vignette("textreuse-minhash", package = "textreuse")
#  vignette("textreuse-alignment", package = "textreuse")

## ------------------------------------------------------------------------
library(textreuse)
file <- system.file("extdata/ats/remember00palm.txt", 
                    package = "textreuse")
doc <- TextReuseTextDocument(file = file, meta = list("publisher" = "ATS"),
                             tokenizer = tokenize_ngrams, n = 5,
                             keep_tokens = TRUE)
doc

## ------------------------------------------------------------------------
meta(doc)
meta(doc, "id")
meta(doc, "date") <- 1865
head(tokens(doc))
head(hashes(doc))
wordcount(doc)

## ------------------------------------------------------------------------
dir <- system.file("extdata/ats", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir, tokenizer = tokenize_ngrams, n = 5,
                          progress = FALSE)
corpus

## ------------------------------------------------------------------------
names(corpus)
corpus[["remember00palm"]]
corpus[c("calltounconv00baxt", "lifeofrevrichard00baxt")]

## ------------------------------------------------------------------------
wordcount(corpus)

## ------------------------------------------------------------------------
text <- "How many roads must a man walk down\nBefore you'll call him a man?"

tokenize_words(text)
tokenize_sentences(text)
tokenize_ngrams(text, n = 3)
tokenize_skip_ngrams(text, n = 3, k = 2)

## ------------------------------------------------------------------------
poem <- "Roses are red\nViolets are blue\nI like using R\nAnd you should too"
cat(poem)

tokenize_lines <- function(string) {
  stringr::str_split(string, "\n+")[[1]]
}

tokenize_lines(poem)

## ------------------------------------------------------------------------
hash_string(tokenize_words(text))

## ------------------------------------------------------------------------
a <- tokenize_words(paste("How does it feel, how does it feel?",
                          "To be without a home",
                          "Like a complete unknown, like a rolling stone"))
b <- tokenize_words(paste("How does it feel, how does it feel?",
                          "To be on your own, with no direction home",
                          "A complete unknown, like a rolling stone"))

jaccard_similarity(a, b)
jaccard_dissimilarity(a, b)
jaccard_bag_similarity(a, b)
ratio_of_matches(a, b)

## ----eval = FALSE--------------------------------------------------------
#  options("mc.cores" = 4L)

