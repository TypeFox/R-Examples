######################
### proximity database
######################

### check functions
.function_or_character <-
function(x)
{
    if (!is.character(x) && !is.function(x))
        stop("Need function or function name.")
}

.abcd_and_binary <-
function(x)
{
    if (x$abcd && (x$type != "binary"))
        stop(paste(dQuote("abcd"), "mode only available for binary measures."))
}

### create registry
pr_DB <- registry(registry_class = "pr_DB",
                  entry_class = "proxy_registry_entry",
                  validity_FUN = ".abcd_and_binary")

## create fields
pr_DB$set_field("FUN", is_mandatory = TRUE, validity_FUN = ".function_or_character")
pr_DB$set_field("distance", type = "logical", default = TRUE)
pr_DB$set_field("PREFUN", validity_FUN = ".function_or_character")
pr_DB$set_field("POSTFUN", validity_FUN = ".function_or_character")
pr_DB$set_field("convert", validity_FUN = ".function_or_character")
pr_DB$set_field("type", type = c("binary", "nominal", "ordinal", "metric", "other"),
                default = "other")
pr_DB$set_field("loop", type = "logical", default = TRUE)
pr_DB$set_field("C_FUN", type = "logical", default = FALSE)
pr_DB$set_field("PACKAGE", type = "character", default = "proxy")
pr_DB$set_field("abcd", type = "logical", default = FALSE)
pr_DB$set_field("formula", type = "character")
pr_DB$set_field("reference", type = "character")
pr_DB$set_field("description", type = "character")

### summary and print methods
summary.pr_DB <-
function(object, verbosity = c("short", "long"), ...)
{
    if (length(object) < 1)
        return(object)
    verbosity <- match.arg(verbosity)

    object <-
        switch(verbosity,
               short = list(names = object$get_field_entries("names", unlist = FALSE),
                            distance = object$get_field_entries("distance")),
               long = list(names = object$get_field_entries("names", unlist = FALSE),
                           distance = object$get_field_entries("distance"),
                           type = object$get_field_entries("type"),
                           formula = object$get_field_entries("formula"))
               )

    structure(object, class = "summary.pr_DB")
}

print.summary.pr_DB <- function(x, ...)
{
    distance <- c("Similarity", "Distance")[x[[2]] + 1]
    if (length(x) > 2)
        x[[3]][is.na(x[[3]])] <- "other"
    for (i in unique(distance)) {
        ind <- which(distance == i)
        if (length(ind) > 0) {
            writeLines(paste("*", i, "measures:"))
            if (length(x) > 2) {
                for (k in ind)
                    writeLines(paste("     ", paste(x[[1]][[k]], collapse = "/"),
                                     " (", x[[3]][k], ") = ", x[[4]][k], sep = ""))
            } else {
                tmp <- sort(sapply(x[[1]][ind], function(i) i[1]))
                writeLines(strwrap(paste(tmp, collapse = ", ")))
            }
        }
        writeLines("")
    }
}
