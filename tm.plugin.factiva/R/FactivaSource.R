FactivaSource <- function(x, encoding = "UTF-8", format = c("auto", "XML", "HTML")) {
    format <- match.arg(format)

    # XML format
    if(format == "XML" ||
       (format == "auto" && grepl(".(xml|XML)$", x))) {
        XMLSource(x, function(tree) {
            sapply(xmlChildren(xmlChildren(xmlChildren(xmlRoot(tree))
                $ppsArticleResponse)$ppsarticleResultSet), xmlChildren)
        },
                      readFactivaXML)
    }
    # HTML format
    else {
        tree <- htmlParse(x, encoding=encoding)

        # The full class is "article XXArticle", with XX the language code
        content <- getNodeSet(tree, "//div[starts-with(@class, 'article')]")
        free(tree)

        SimpleSource(encoding, length=length(content),
                     content=content, uri=x,
                     reader=readFactivaHTML, class="FactivaSource")
    }
}

# This functions need to be exactly the same as those for XMLSource
# since they can be used with the Factiva XML source
getElem.FactivaSource <- function(x) list(content = saveXML(x$content[[x$position]]), uri = x$uri)
