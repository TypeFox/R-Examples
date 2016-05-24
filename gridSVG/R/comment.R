grid.comment <- function(comment, name = NULL, vp = NULL) {
    grid.draw(commentGrob(comment, name, vp))
}

commentGrob <- function(comment, name = NULL, vp = NULL) {
    g <- grob(name = name, vp = vp, cl = "comment")
    g$comment <- comment
    cl <- class(g)
    class(g) <- unique(c("comment.grob", cl))
    g
}

primToDev.comment.grob <- function(x, dev) {
    svgComment(x$comment, dev@dev)
}

