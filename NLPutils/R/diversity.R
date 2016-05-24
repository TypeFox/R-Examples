QDAP_Diversity_Annotator <-
function()
{
    f <- function(s, a) {
        s <- as.String(s)
        diversity <- tail(unclass(diversity(s)), -2L)
        Annotation(.document_id(a),
                   "document",
                   1L,
                   nchar(s),
                   list(list(diversity = diversity)))
    }

    description <-
        "Document diversity annotator for English using qdap::diversity"

    Annotator(f, list(description = description))
}
