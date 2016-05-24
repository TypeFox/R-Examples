# Functions for generating arbitrary SVG elements

elementGrob <- function(el, name = NULL, attrs = NULL, 
                        namespace = NULL,
                        namespaceDefinitions = NULL,
                        children = NULL,
                        vp = NULL, childrenvp = NULL,
                        asis = FALSE) {
    eg <- gTree(name = name, vp = vp,
                children = children, childrenvp = childrenvp,
                cl = "element")
    # Keeping copy of name because of asis.
    # If it's TRUE, we leave the id alone.
    # When FALSE, the resulting id attribute could get modified
    # by things like gTrees so that the name is a *path*.
    eg$asis <- asis
    eg$origname <- eg$name

    eg$el <- el
    eg$attrs <- if (is.null(attrs)) list() else attrs
    eg$namespace <- namespace
    eg$namespaceDefinitions <- namespaceDefinitions
    cl <- class(eg)
    class(eg) <- unique(c("element.grob", cl))
    eg
}

grid.element <- function(el, name = NULL, attrs = NULL, 
                         namespace = NULL,
                         namespaceDefinitions = NULL,
                         children = NULL,
                         vp = NULL, childrenvp = NULL,
                         asis = FALSE) {
    grid.draw(elementGrob(el, name, attrs, namespace,
                          namespaceDefinitions, children,
                          vp, childrenvp, asis))
}

devGrob.element.grob <- function(x, dev) {
  list(id = if (x$asis) x$origname
            else getID(x$name, "grob"),
       name = x$el,
       classes = x$classes,
       attrs = x$attrs,
       namespace = x$namespace,
       namespaceDefinitions = x$namespaceDefinitions)
}

# Unlike gTrees, we don't need a group for children because it
# complicates output, when we want clear output to SVG.
# Also, do *not* add gpars because they also complicate output,
# if we *really* want to do it, then just use the 'attrs' arg.
primToDev.element.grob <- function(x, dev) {
    devStartElement(devGrob(x, dev), NULL, dev)
    lapply(x$children, function(child) {
        grobToDev(child, dev)
    })
    devEndElement(x$name, dev)
}

devGrob.textnode.grob <- function(x, dev) {
    list(text = x$text)
}

primToDev.textnode.grob <- function(x, dev) {
    devTextNode(devGrob(x, dev), dev)
}

textNodeGrob <- function(text, name = NULL, vp = NULL) {
    if (length(text) > 1)
        stop("'text' must be a single element character vector")
    tng <- grob(name = name, vp = vp, cl = "textNode")
    tng$text <- text
    cl <- class(tng)
    class(tng) <- unique(c("textnode.grob", cl))
    tng
}

grid.textNode <- function(text, name = NULL, vp = NULL) {
    grid.draw(textNodeGrob(text, name, vp))
}
