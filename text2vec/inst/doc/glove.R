## ---- eval=FALSE---------------------------------------------------------
#  library(text2vec)
#  library(readr)
#  temp <- tempfile()
#  download.file('http://mattmahoney.net/dc/text8.zip', temp)
#  wiki <- read_lines(unz(temp, "text8"))
#  unlink(temp)

## ---- eval=FALSE---------------------------------------------------------
#  # Create iterator over tokens
#  tokens <- strsplit(wiki, split = " ", fixed = T)
#  # Create vocabulary. Terms will be unigrams (simple words).
#  vocab <- create_vocabulary(itoken(tokens))

## ---- eval=FALSE---------------------------------------------------------
#  vocab <- prune_vocabulary(vocab, term_count_min = 5L)

## ---- eval=FALSE---------------------------------------------------------
#  # We provide an iterator to create_vocab_corpus function
#  it <- itoken(tokens)
#  # Use our filtered vocabulary
#  vectorizer <- vocab_vectorizer(vocab,
#                                 # don't vectorize input
#                                 grow_dtm = FALSE,
#                                 # use window of 5 for context words
#                                 skip_grams_window = 5L)
#  tcm <- create_tcm(it, vectorizer)

## ---- eval = FALSE-------------------------------------------------------
#  fit <- glove(tcm = tcm,
#               word_vectors_size = 50,
#               x_max = 10, learning_rate = 0.2,
#               num_iters = 15)

## ---- eval = FALSE-------------------------------------------------------
#  word_vectors <- fit$word_vectors[[1]] + fit$word_vectors[[2]]
#  rownames(word_vectors) <- rownames(tcm)

## ---- eval = FALSE-------------------------------------------------------
#  word_vectors_norm <- sqrt(rowSums(word_vectors ^ 2))
#  
#  rome <- word_vectors['paris', , drop = FALSE] -
#    word_vectors['france', , drop = FALSE] +
#    word_vectors['italy', , drop = FALSE]
#  
#  cos_dist <- text2vec:::cosine(rome,
#                                word_vectors,
#                                word_vectors_norm)
#  head(sort(cos_dist[1,], decreasing = T), 10)
#  ##     paris    venice     genoa      rome  florence
#  ## 0.7811252 0.7763088 0.7048109 0.6696540 0.6580989

