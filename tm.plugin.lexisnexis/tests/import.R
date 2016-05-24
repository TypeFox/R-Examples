library(tm)
library(tm.plugin.lexisnexis)

# English example
file <- system.file("texts", "lexisnexis_test_en.html",
                    package = "tm.plugin.lexisnexis")
corpus <- Corpus(LexisNexisSource(file))

# French example
file <- system.file("texts", "lexisnexis_test_fr.html",
                    package = "tm.plugin.lexisnexis")
corpus <- Corpus(LexisNexisSource(file))
