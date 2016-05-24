emotions <- mldr_from_dataframe(mldr::emotions$dataset[, c(emotions$attributesIndexes, emotions$labels$index)], mldr::emotions$labels$index, "emotions")
birds <- mldr_from_dataframe(mldr::birds$dataset[, c(birds$attributesIndexes, birds$labels$index)], mldr::birds$labels$index, "birds")
genbase <- mldr_from_dataframe(mldr::genbase$dataset[, c(genbase$attributesIndexes, genbase$labels$index)], mldr::genbase$labels$index, "genbase")

save(emotions, file = "emotions.rda", compress = "xz")
save(birds, file = "birds.rda", compress = "xz")
save(genbase, file = "genbase.rda", compress = "xz")
