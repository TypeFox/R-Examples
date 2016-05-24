QDAP_Formality_Annotator <-
function()
{
    f <- function(s, a) {
        s <- as.String(s)
        formality <- formality(s, progress.bar = FALSE)$formality$formality
        Annotation(.document_id(a),
                   "document",
                   1L,
                   nchar(s),
                   list(list(formality = formality)))
    }

    description <-
        "Document formality annotator for English using qdap::formality"

    Annotator(f, list(description = description))
}
