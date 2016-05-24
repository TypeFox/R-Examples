MalletLDA <- function(num.topics = 10, alpha.sum = 5.0, beta = 0.01) { .jnew("cc/mallet/topics/RTopicModel", num.topics, alpha.sum, beta) }

mallet.topic.words <- function(topic.model, normalized=FALSE, smoothed=FALSE) { .jevalArray(topic.model$getTopicWords(normalized, smoothed), simplify=T) }
mallet.doc.topics <- function(topic.model, normalized=FALSE, smoothed=FALSE) { .jevalArray(topic.model$getDocumentTopics(normalized, smoothed), simplify=T) }

mallet.word.freqs <- function(topic.model) {
  word.freqs <- .jevalArray(topic.model$getWordFrequencies(), simplify=T)
  data.frame(words = topic.model$getVocabulary(), term.freq = word.freqs[,1], doc.freq = word.freqs[,2])
}

mallet.subset.topic.words <- function(topic.model, subset.docs, normalized=FALSE, smoothed=FALSE) {
  .jevalArray(topic.model$getSubCorpusTopicWords(subset.docs, normalized, smoothed), simplify=T)
}

mallet.top.words <- function(topic.model, word.weights, num.top.words=10) {
  top.indices <- order(word.weights, decreasing=T)[1:num.top.words]
  data.frame(words = topic.model$getVocabulary()[top.indices], weights = word.weights[top.indices], stringsAsFactors=F)
}

mallet.import <- function(id.array, text.array, stoplist.file, preserve.case=FALSE, token.regexp="[\\p{L}]+") {
  stoplist.file <- normalizePath(stoplist.file)
  if (class(text.array[1]) != "character") stop("Text field is not a string. Remember to create data frames with stringsAsFactors=F.")
  token.pattern <- J("java/util/regex/Pattern")$compile(token.regexp)
  pipe.list <- .jnew("java/util/ArrayList")
  pipe.list$add(.jnew("cc/mallet/pipe/CharSequence2TokenSequence", token.pattern))
  if (! preserve.case) { pipe.list$add(.jnew("cc/mallet/pipe/TokenSequenceLowercase")) }
  pipe.list$add(.jnew("cc/mallet/pipe/TokenSequenceRemoveStopwords", .jnew("java/io/File", stoplist.file), "UTF-8", FALSE, FALSE, FALSE))
  pipe.list$add(.jnew("cc/mallet/pipe/TokenSequence2FeatureSequence"))
  #pipe.list$add(.jnew("cc/mallet/pipe/PrintInputAndTarget"))

  pipe <- .jnew("cc/mallet/pipe/SerialPipes", .jcast(pipe.list, "java/util/Collection"))

  instances <- .jnew("cc/mallet/types/InstanceList", .jcast(pipe, "cc/mallet/pipe/Pipe"))

  J("cc/mallet/topics/RTopicModel")$addInstances(instances, id.array, text.array)

  return(instances)
}

# mallet.read.dir() function, created by Dan Bowen
# This function takes a directory path as its only argument
# ... and returns a data.frame() with 2 columns: <id> & <text>.
# ... This data.frame() has as many rows as there are files in the Dir.
# The form of this functions return attempts to conform to that
# ... used by the mallet.import() function, available in the 'mallet' R package
mallet.read.dir <- function(Dir) {
  # get Dir Files (filepaths)
  Files <- file.path(Dir, list.files(Dir))
  # for each File:
  mallet.read.file <- function(File) {
    # read File, per line
    Lines <- scan(File, what='character', sep='\n', quote='')
    # paste Lines back together with '\n'
    string <- paste(Lines, collapse='\n')
    # return data.frame
    data.frame(id=File, text=string, stringsAsFactors=F)
  }
  # apply the above function to the Files in the dir
  # ... rbind the resulting list of data.frames together
  do.call(rbind, lapply(Files, mallet.read.file))
}

## Get a vector containing short names for all the topics
mallet.topic.labels <- function(topic.model, topic.words, num.top.words=3) {
  n.topics <- dim(topic.words)[1]
  topics.labels <- rep("", n.topics)
  for (topic in 1:n.topics) topics.labels[topic] <- paste(mallet.top.words(topic.model, topic.words[topic,], num.top.words=3)$words, collapse=" ")
  topics.labels
}

## Return a hierarchical clustering of topics.
mallet.topic.hclust <- function(doc.topics, topic.words, balance = 0.3) {
  ## transpose and normalize the doc topics
  topic.docs <- t(doc.topics)
  topic.docs <- topic.docs / rowSums(topic.docs)

  hclust(balance * dist(topic.words) + (1.0 - balance) * dist(topic.docs))
}