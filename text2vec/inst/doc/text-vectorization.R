## ---- echo=FALSE--------------------------------------------------------------
op <- options(width = 80, str = strOptions(strict.width = "cut"))

## ---- loading-data, eval=TRUE-------------------------------------------------
library(text2vec)
data("movie_review")
set.seed(42L)

## ---- vocab-iterator, eval=TRUE-----------------------------------------------
it <- itoken(movie_review$review, 
             preprocess_function = tolower, 
             tokenizer = word_tokenizer, 
             ids = movie_review$id)

sw <- c("i", "me", "my", "myself", "we", "our", "ours", "ourselves", "you", "your", "yours")
vocab <- create_vocabulary(it, stopwords = sw)

## -----------------------------------------------------------------------------
# Each element of list represents document
tokens <- movie_review$review %>% 
  tolower() %>% 
  word_tokenizer()
it <- itoken(tokens, ids = movie_review$id)
vocab <- create_vocabulary(it, stopwords = sw)

## ---- vocab_dtm_1, eval=TRUE--------------------------------------------------
it <- itoken(tokens, ids = movie_review$id)
# Or
# it <- itoken(movie_review$review, tolower, word_tokenizer, ids = movie_review$id)
vectorizer <- vocab_vectorizer(vocab)
dtm <- create_dtm(it, vectorizer)

## ---- vocab_dtm_1_dim, eval=TRUE----------------------------------------------
str(dtm)
identical(rownames(dtm), movie_review$id)

## ---- fit_1, message=FALSE, warning=FALSE, eval=TRUE--------------------------
library(glmnet)
fit <- cv.glmnet(x = dtm, y = movie_review[['sentiment']], 
                 family = 'binomial', 
                 # lasso penalty
                 alpha = 1,
                 # interested in the area under ROC curve
                 type.measure = "auc",
                 # 5-fold cross-validation
                 nfolds = 5,
                 # high value is less accurate, but has faster training
                 thresh = 1e-3,
                 # again lower number of iterations for faster training
                 maxit = 1e3)
plot(fit)
print(paste("max AUC =", round(max(fit$cvm), 4)))

## ---- prune_vocab_dtm_1-------------------------------------------------------
pruned_vocab <- prune_vocabulary(vocab, term_count_min = 10,
 doc_proportion_max = 0.5, doc_proportion_min = 0.001)
it <- itoken(tokens, ids = movie_review$id)
vectorizer <- vocab_vectorizer(pruned_vocab)
dtm <- create_dtm(it, vectorizer)
dim(dtm)

## ---- tfidf_dtm_1-------------------------------------------------------------
dtm <- dtm %>% transform_tfidf()

## ---- fit_2, message=FALSE, warning=FALSE, eval=TRUE--------------------------
t1 <- Sys.time()
fit <- cv.glmnet(x = dtm, y = movie_review[['sentiment']], 
                 family = 'binomial', 
                 alpha = 1,
                 type.measure = "auc",
                 nfolds = 5,
                 thresh = 1e-3,
                 maxit = 1e3)
print(difftime(Sys.time(), t1, units = 'sec'))
plot(fit)
print(paste("max AUC =", round(max(fit$cvm), 4)))

## ---- ngram_dtm_1-------------------------------------------------------------
it <- itoken(tokens, ids = movie_review$id)

vocab <- create_vocabulary(it, ngram = c(1L, 3L)) %>% 
  prune_vocabulary(term_count_min = 10, 
                   doc_proportion_max = 0.5, 
                   doc_proportion_min = 0.001)

vectorizer <- vocab_vectorizer(vocab)

dtm <- tokens %>% 
  itoken() %>% 
  create_dtm(vectorizer) %>% 
  transform_tfidf()

dim(dtm)

fit <- cv.glmnet(x = dtm, y = movie_review[['sentiment']], 
                 family = 'binomial', 
                 alpha = 1,
                 type.measure = "auc",
                 nfolds = 5,
                 thresh = 1e-3,
                 maxit = 1e3)

plot(fit)
print(paste("max AUC =", round(max(fit$cvm), 4)))

## ---- hash_dtm----------------------------------------------------------------
it <- itoken(tokens, ids = movie_review$id)

vectorizer <- hash_vectorizer(hash_size = 2 ^ 16, ngram = c(1L, 3L))
dtm <- create_dtm(it, vectorizer) %>% 
  transform_tfidf()

fit <- cv.glmnet(x = dtm, y = movie_review[['sentiment']], 
                 family = 'binomial', 
                 alpha = 1,
                 type.measure = "auc",
                 nfolds = 5,
                 thresh = 1e-3,
                 maxit = 1e3)

plot(fit)
print(paste("max AUC =", round(max(fit$cvm), 4)))

