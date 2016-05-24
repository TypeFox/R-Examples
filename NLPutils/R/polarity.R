QDAP_Polarity_Annotator <-
function()
{
    f <- function(s, a) {
        s <- as.String(s)
        a <- a[a$type == "sentence"]
        p <- polarity(s[a])$all$polarity
        a$features <- lapply(p, single_feature, "polarity")
        a
    }

    description <-
        "Sentence polarity annotator for English using qdap::polarity"

    Annotator(f, list(description = description))
}
