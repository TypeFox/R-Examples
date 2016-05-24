# Bouchet-Valat & Bastin, RcmdrPlugin.temis, a Graphical Integrated Text Mining Solution in R,
# R Journal, 5(1), 2013.

# Stopwords change accross tm versions
stopwords.en <- c("a", "about", "above", "across", "after", "again", 
"against", "all", "almost", "alone", "along", "already", "also", 
"although", "always", "am", "among", "an", "and", "another", 
"any", "anybody", "anyone", "anything", "anywhere", "are", "area", 
"areas", "aren't", "around", "as", "ask", "asked", "asking", 
"asks", "at", "away", "b", "back", "backed", "backing", "backs", 
"be", "became", "because", "become", "becomes", "been", "before", 
"began", "behind", "being", "beings", "below", "best", "better", 
"between", "big", "both", "but", "by", "c", "came", "can", "can't", 
"cannot", "case", "cases", "certain", "certainly", "clear", "clearly", 
"come", "could", "couldn't", "d", "did", "didn't", "differ", 
"different", "differently", "do", "does", "doesn't", "doing", 
"don't", "done", "down", "downed", "downing", "downs", "during", 
"e", "each", "early", "either", "end", "ended", "ending", "ends", 
"enough", "even", "evenly", "ever", "every", "everybody", "everyone", 
"everything", "everywhere", "f", "face", "faces", "fact", "facts", 
"far", "felt", "few", "find", "finds", "first", "for", "four", 
"from", "full", "fully", "further", "furthered", "furthering", 
"furthers", "g", "gave", "general", "generally", "get", "gets", 
"give", "given", "gives", "go", "going", "good", "goods", "got", 
"great", "greater", "greatest", "group", "grouped", "grouping", 
"groups", "h", "had", "hadn't", "has", "hasn't", "have", "haven't", 
"having", "he", "he'd", "he'll", "he's", "her", "here", "here's", 
"hers", "herself", "high", "higher", "highest", "him", "himself", 
"his", "how", "how's", "however", "i", "i'd", "i'll", "i'm", 
"i've", "if", "important", "in", "interest", "interested", "interesting", 
"interests", "into", "is", "isn't", "it", "it's", "its", "itself", 
"j", "just", "k", "keep", "keeps", "kind", "knew", "know", "known", 
"knows", "l", "large", "largely", "last", "later", "latest", 
"least", "less", "let", "let's", "lets", "like", "likely", "long", 
"longer", "longest", "m", "made", "make", "making", "man", "many", 
"may", "me", "member", "members", "men", "might", "more", "most", 
"mostly", "mr", "mrs", "much", "must", "mustn't", "my", "myself", 
"n", "necessary", "need", "needed", "needing", "needs", "never", 
"new", "newer", "newest", "next", "no", "nobody", "non", "noone", 
"nor", "not", "nothing", "now", "nowhere", "number", "numbers", 
"o", "of", "off", "often", "old", "older", "oldest", "on", "once", 
"one", "only", "open", "opened", "opening", "opens", "or", "order", 
"ordered", "ordering", "orders", "other", "others", "ought", 
"our", "ours", "ourselves", "out", "over", "own", "p", "part", 
"parted", "parting", "parts", "per", "perhaps", "place", "places", 
"point", "pointed", "pointing", "points", "possible", "present", 
"presented", "presenting", "presents", "problem", "problems", 
"put", "puts", "q", "quite", "r", "rather", "really", "right", 
"room", "rooms", "s", "said", "same", "saw", "say", "says", "second", 
"seconds", "see", "seem", "seemed", "seeming", "seems", "sees", 
"several", "shall", "shan't", "she", "she'd", "she'll", "she's", 
"should", "shouldn't", "show", "showed", "showing", "shows", 
"side", "sides", "since", "small", "smaller", "smallest", "so", 
"some", "somebody", "someone", "something", "somewhere", "state", 
"states", "still", "such", "sure", "t", "take", "taken", "than", 
"that", "that's", "the", "their", "theirs", "them", "themselves", 
"then", "there", "there's", "therefore", "these", "they", "they'd", 
"they'll", "they're", "they've", "thing", "things", "think", 
"thinks", "this", "those", "though", "thought", "thoughts", "three", 
"through", "thus", "to", "today", "together", "too", "took", 
"toward", "turn", "turned", "turning", "turns", "two", "u", "under", 
"until", "up", "upon", "us", "use", "used", "uses", "v", "very", 
"w", "want", "wanted", "wanting", "wants", "was", "wasn't", "way", 
"ways", "we", "we'd", "we'll", "we're", "we've", "well", "wells", 
"went", "were", "weren't", "what", "what's", "when", "when's", 
"where", "where's", "whether", "which", "while", "who", "who's", 
"whole", "whom", "whose", "why", "why's", "will", "with", "within", 
"without", "won't", "work", "worked", "working", "works", "would", 
"wouldn't", "x", "y", "year", "years", "yes", "yet", "you", "you'd", 
"you'll", "you're", "you've", "young", "younger", "youngest", 
"your", "yours", "yourself", "yourselves", "z")

library(RcmdrPlugin.temis)
library(tm.plugin.factiva)
library(SnowballC)

corpus <- Corpus(FactivaSource(system.file("texts", "reut21578-factiva.xml", package="tm.plugin.factiva")), 
                 readerControl=list(language="en"))
names(corpus) <- make.unique(names(corpus))
corpusVars <- extractMetadata(corpus)
corpusVars <- corpusVars[c("Origin", "Date", "United.States", "North.America", "Canada", 
  "Ecuador", "South.America", "Kuwait", "Middle.East", "Indonesia", "Asia", "Bahrain", 
  "Saudi.Arabia", "Qatar", "United.Arab.Emirates", "Argentina")]
meta(corpus, "Date") <- corpusVars$Date
dtmCorpus <- corpus
dtmCorpus <- tm_map(dtmCorpus, content_transformer(tolower))
dtmCorpus <- tm_map(dtmCorpus, content_transformer(function(x) 
  gsub("(['<U+2019>\n<U+202F><U+2009>]|[[:punct:]]|[[:space:]]|[[:cntrl:]])+", " ", x)))
dtmCorpus <- tm_map(dtmCorpus, removeNumbers)
dtm <- DocumentTermMatrix(dtmCorpus, control=list(tolower=FALSE, wordLengths=c(2, Inf)))
rm(dtmCorpus)
dictionary <- data.frame(row.names=colnames(dtm), "Occurrences"=col_sums(dtm), 
  "Stemmed.Term"=wordStem(colnames(dtm), "en"), "Stopword"=ifelse(colnames(dtm) %in% 
  stopwords("en"), "Stopword", ""), stringsAsFactors=FALSE)
dtm <- dtm[, !colnames(dtm) %in% stopwords.en]
dtm <- rollup(dtm, 2, dictionary[colnames(dtm), 2])
attr(dtm, "dictionary") <- dictionary
rm(dictionary)
meta(corpus, type="corpus", tag="language") <- attr(dtm, "language") <- "en"
meta(corpus, type="corpus", tag="processing") <- attr(dtm, "processing") <- 
  c(lowercase=TRUE, punctuation=TRUE, digits=TRUE, stopwords=TRUE, stemming=TRUE, 
  customStemming=FALSE, twitter=FALSE, removeHashtags=NA, removeNames=NA)
corpus
dtm

# Table 2
specTerms <- specificTerms(dtm, meta(corpus, "Date")[[1]], p=0.1, min.occ=5, n.max=25)
attr(specTerms, "title") <- "Specific terms by Date"
stopifnot(all.equal(round(c(specTerms[[1]]), 4),
c(2.0101, 1.2563, 1.5075, 1.7588, 1.005, 1.2563, 2.5126, 
1.7588, 4.0201, 1.005, 0.7538, 0.7538, 0.7538, 1.005, 1.005, 
4.2714, 5.5276, 0.7538, 0.7538, 0, 0, 0, 0, 0, 0, 66.6667, 71.4286, 
60, 50, 66.6667, 55.5556, 38.4615, 43.75, 30.7692, 50, 60, 60, 
60, 44.4444, 44.4444, 26.9841, 25.5814, 50, 50, 0, 0, 0, 0, 0, 
0, 0.5794, 0.338, 0.4829, 0.676, 0.2897, 0.4346, 1.2554, 0.7726, 
2.5109, 0.3863, 0.2414, 0.2414, 0.2414, 0.4346, 0.4346, 3.042, 
4.1526, 0.2897, 0.2897, 1.0623, 0.9174, 0.8209, 0.676, 0.6277, 
0.5311, 8, 5, 6, 7, 4, 5, 10, 7, 16, 4, 3, 3, 3, 4, 4, 17, 22, 
3, 3, 0, 0, 0, 0, 0, 0, 12, 7, 10, 14, 6, 9, 26, 16, 52, 8, 5, 
5, 5, 9, 9, 63, 86, 6, 6, 22, 19, 17, 14, 13, 11, 3.3408, 2.6673, 
2.5714, 2.3631, 2.1825, 2.1365, 2.1162, 2.0293, 1.8867, 1.6542, 
1.627, 1.627, 1.627, 1.4379, 1.4379, 1.4006, 1.3699, 1.3456, 
1.3456, -2.3698, -2.1201, -1.9406, -1.6461, -1.5396, -1.3106, 
0.0004, 0.0038, 0.0051, 0.0091, 0.0145, 0.0163, 0.0172, 0.0212, 
0.0296, 0.049, 0.0519, 0.0519, 0.0519, 0.0752, 0.0752, 0.0807, 
0.0854, 0.0892, 0.0892, 0.0089, 0.017, 0.0262, 0.0499, 0.0618, 
0.095), check.attributes=FALSE))

# Table 3
dissDtm <- rollup(dtm, 1, meta(corpus, "Date"))
diss <- dist(sweep(dissDtm/row_sums(dissDtm), 2, sqrt(sum(dissDtm)/col_sums(dissDtm)), 
  "*"))
rm(dissDtm)
attr(diss, "title") <- "Date by Date dissimilarity table"
diss

