handlersSDML <- function() {
    list(
         textdata = function(x, ...) {
             
             sep <- ifelse("sep" %in% names(xmlAttrs(x)),
                           xmlAttrs(x)["sep"], " \n\r")
             sep <- paste("[", sep, "]+", sep="")
             
             type <- ifelse("type" %in% names(xmlAttrs(x)),
                            xmlAttrs(x)["type"], "character")
             mode <- NULL
             if ("mode" %in% names(xmlAttrs(x)))
                 mode <- as.character(xmlAttrs(x)["mode"])
             posinf <- ifelse("posinf.string" %in% names(xmlAttrs(x)),
                              xmlAttrs(x)["posinf.string"], "+Inf")
             neginf <- ifelse("neginf.string" %in% names(xmlAttrs(x)),
                              xmlAttrs(x)["neginf.string"],"-Inf")
             nan <- ifelse("nan.string" %in% names(xmlAttrs(x)),
                           xmlAttrs(x)["nan.string"], "NaN")
             na <- ifelse("na.string" %in% names(xmlAttrs(x)),
                          xmlAttrs(x)["na.string"], "NA")
             null <- ifelse("null.string" %in% names(xmlAttrs(x)),
                            xmlAttrs(x)["null.string"], "NULL")
             true <- ifelse("true" %in% names(xmlAttrs(x)),
                            xmlAttrs(x)["true"], "1")
             false <- ifelse("false" %in% names(xmlAttrs(x)),
                             xmlAttrs(x)["false"], "0")
             
             children <- if(is.null(x[[1]]))
                 NULL
             else {
                 tmp <- unlist(strsplit(xmlValue(x[[1]]), sep))
                 list(tmp[tmp != ""])
             }
             
             structure(list(children = children,
                            name = "textdata",
                            attributes = c(
                            type = type,
                            mode = mode,
                            posinf.string = posinf,
                            neginf.string = neginf,
                            nan.string = nan,
                            na.string = na,
                            null.string = null,
                            true = true,
                            false = false
                            )
                            ),
                       class = "XMLNode"
                       )
         },
         e = function(x, ...) {
             info <- ifelse("info" %in% names(xmlAttrs(x)),
                            xmlAttrs(x)["info"], NA)
             children <- if(!length(xmlChildren(x)))
                 NULL
             else
                 xmlChildren(x)
             structure(list(children = children,
                            name = "e",
                            attributes = c(
                            info = info
                            )
                            ),
                       class = "XMLNode"
                       )
         }
         )
}

