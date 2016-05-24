read_tag_dictionary <-
function(dict) {
    it <- dict$iterator()
    tags <- words <- list()
    while(it$hasNext()) {
        w <- .jcall(it, "Ljava/lang/Object;", "next")
        words <- c(words, list(.jsimplify(w)))
        tags <- c(tags, list(dict$getTags(w)))
    }
    names(tags) <- unlist(words)
    tags
}
