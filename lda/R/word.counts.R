word.counts <-
function (docs, vocab = NULL) 
{
    all.words <- matrix(unlist(docs), nrow = 2)
    if (is.null(vocab)) {
        tapply(all.words[2, ], all.words[1, ], sum)
    }
    else {
        result <- tapply(all.words[2, ], list(word = ordered(vocab[all.words[1, 
            ] + 1], levels = vocab)), sum)
        result[is.na(result)] <- 0
        result
    }
}
