## ------------------------------------------------------------------------
library(text2vec)
library(magrittr)
data("movie_review")
# remove all internal EOL to simplify reading
movie_review$review <- gsub(pattern = '\n', replacement = ' ', 
                            x = movie_review$review, fixed = TRUE)
N_FILES <- 10
CHUNK_LEN <- nrow(movie_review) / N_FILES
files <- sapply(1:N_FILES, function(x) tempfile())
chunks <- split(movie_review, rep(1:N_FILES, each = nrow(movie_review) / N_FILES ))
for (i in 1:N_FILES ) {
  write.table(chunks[[i]], files[[i]], quote = T, row.names = F, col.names = T, sep = '|')
}
# Note how data looks like
str(movie_review, strict.width = 'cut')

## ------------------------------------------------------------------------
library(data.table)
reader <- function(x, ...) {
  # read
  chunk <- fread(x, header = T, sep = '|')
  # select column with review
  res <- chunk$review
  # assign ids to reviews
  names(res) <- chunk$id
  res
}
# create iterator over files
it_files  <- ifiles(files, reader_function = reader)
# create iterator over tokens from files iterator
it_tokens <- itoken(it_files, preprocess_function = tolower, 
                    tokenizer = word_tokenizer, progessbar = FALSE)

vocab <- create_vocabulary(it_tokens)

## ------------------------------------------------------------------------
# need to reinitialise iterators!
# they are mutable and already empty!
# try(it_files$nextElem())
it_files  <- ifiles(files, reader_function = reader)
it_tokens <- itoken(it_files, preprocess_function = tolower, 
                    tokenizer = word_tokenizer, progessbar = FALSE)

dtm <- create_dtm(it_tokens, vectorizer = vocab_vectorizer(vocab), type = 'lda_c')
str(dtm, list.len = 5)

## ---- eval=FALSE---------------------------------------------------------
#  library(lda)
#  # prior for topics
#  alpha = 0.1
#  # prior for words
#  eta = 0.001
#  # fit model with 30 topics, make 30 Gibbs sampling iterations
#  lda_fit <- lda.collapsed.gibbs.sampler(documents = dtm, K = 30,
#                                         vocab = vocab$vocab$terms,
#                                         alpha = alpha,
#                                         eta = eta,
#                                         num.iterations = 30,
#                                         trace = 2L)

## ---- warning=FALSE, message=FALSE, eval=FALSE---------------------------
#  N_WORKERS <- 4
#  library(doParallel)
#  # register parallel backend
#  registerDoParallel(N_WORKERS)
#  
#  #  prepare splits
#  # "jobs" is a list of itoken iterators!
#  N_SPLITS <- 4
#  
#  jobs <- files %>%
#    split_into(N_SPLITS) %>%
#    lapply(ifiles, reader_function = reader) %>%
#    # Worth to set chunks_number to 1 because we already splitted input
#    lapply(itoken, chunks_number = 1, preprocess_function = tolower,
#           tokenizer = word_tokenizer, progessbar = FALSE)
#  
#  # Alternatively when data is in memory we can perform splite in the following way:
#  #
#  # review_chunks <- split_into(movie_review$review, N_SPLITS)
#  # review_ids <- split_into(movie_review$id, N_SPLITS)
#  #
#  # jobs <- Map(function(doc, ids) {
#  #  itoken(iterable = doc, ids = ids, preprocess_function = tolower,
#  #         tokenizer = word_tokenizer, chunks_number = 1, progessbar = FALSE)
#  # }, review_chunks, review_ids)
#  
#  # Now all below function calls will benefit from multicore machines
#  # Each job will be evaluated in separate process
#  
#  # vocabulary creation
#  vocab <- create_vocabulary(jobs)
#  
#  # dtm vocabulary vectorization
#  v_vectorizer <- vocab_vectorizer(vocab)
#  vocab_dtm_parallel <- create_dtm(jobs, vectorizer = v_vectorizer)
#  
#  # dtm hash vectorization
#  h_vectorizer <- hash_vectorizer()
#  hash_dtm_parallel <- create_dtm(jobs, vectorizer = h_vectorizer)
#  
#  # coocurence statistics
#  tcm_vectorizer <- vocab_vectorizer(vocab, grow_dtm = FALSE, skip_grams_window = 5)
#  tcm_parallel <- create_tcm(jobs, vectorizer = tcm_vectorizer)

