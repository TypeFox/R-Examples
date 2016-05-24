library(text2vec)
library(readr)
library(fastmatch)
library(Matrix)
library(stringr)
library(magrittr)


fit <- readRDS("/Users/dmitryselivanov/projects/word_embeddings/enwiki_dim=300_vocab=30k_tcm_pruned/glove_fit.rds")
vocab <- readRDS("/Users/dmitryselivanov/projects/word_embeddings/enwiki_dim=300_vocab=30k/pruned_vocab_text2vec.rds")
terms <- vocab$vocab$terms
word_vectors <- fit$word_vectors[[1]] + fit$word_vectors[[2]]
rownames(word_vectors) <- terms

questions_file_path <- '~/Downloads/datasets/questions-words.txt'
qlst <- prepare_analog_questions(questions_file_path, terms)

Sys.time()
df <- check_accuracy(questions_lst = qlst, m_word_vectors = word_vectors)
Sys.time()

text8 <- read_lines('~/Downloads/datasets/text8')
it <- itoken(text8, preprocess_function = identity, tokenizer = function(x) str_split(x, fixed(" ")))
vocab <- vocabulary(it) %>%
  prune_vocabulary(term_count_min = 5)

it <- itoken(text8, preprocess_function = identity, tokenizer = function(x) str_split(x, fixed(" ")))
vectorizer <- vocab_vectorizer(vocab, grow_dtm = F, skip_grams_window = 5L)
tcm <- create_tcm(it, vectorizer)
# corpus <- create_vocab_corpus(iterator = it, vocabulary = vocab, grow_dtm = F, skip_grams_window = 5)
# tcm <- get_tcm(corpus)

# save(tcm, file = '~/Downloads/datasets/tcm8.RData', compress = F)
# load('~/Downloads/datasets/tcm8.RData')
# load("~/Downloads/m8.RData")
RcppParallel::setThreadOptions(numThreads = 4)
fit <- glove(tcm = tcm, shuffle_seed = 42L, word_vectors_size = 100,  x_max = 10, learning_rate = 0.2,
              num_iters = 10, grain_size = 1e5, max_cost = 10, convergence_threshold = 0.01)

terms <- vocab$vocab$terms
word_vectors <- fit$word_vectors[[1]] + fit$word_vectors[[2]]
rownames(word_vectors) <- terms

questions_file_path <- '~/Downloads/datasets/questions-words.txt'
qlst <- prepare_analogue_questions(questions_file_path, terms)

df <- check_analogue_accuracy(questions_lst = qlst, m_word_vectors = word_vectors)


# top_similar(original_word_wectors['king',], original_word_wectors, 10)
