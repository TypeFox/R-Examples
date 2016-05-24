write.gml <-
function (x, file = NULL) 
{
    if (isTRUE(is.na(dim(x)[3]) == TRUE) == FALSE) {
        if (isTRUE(is.null(dimnames(x)[[3]])) == TRUE) 
            dimnames(x)[[3]] <- 1:dim(x)[3]
    }
    if (isTRUE(is.null(dimnames(x)[[1]])) == TRUE) 
        dimnames(x)[[1]] <- dimnames(x)[[2]] <- 1:dim(x)[1]
    suppressWarnings(file.remove(file = file))
    cat(paste("Creator", "\"multiplex\"", sep = "\t"), file = file, 
        sep = "\n", append = TRUE)
    cat(paste("Version", paste("\"", utils::packageDescription("multiplex")["Version"]$Version, 
        "\"", sep = ""), sep = "\t"), file = file, sep = "\n", 
        append = TRUE)
    cat("graph", file = file, sep = "\n", append = TRUE)
    cat("[", file = file, sep = "\n", append = TRUE)
    cat(paste("", "hierarchic", "1", sep = "\t"), file = file, 
        sep = "\n", append = TRUE)
    cat(paste("", "label", "\"\"", sep = "\t"), file = file, 
        sep = "\n", append = TRUE)
    cat(paste("", "directed", "1", sep = "\t"), file = file, 
        sep = "\n", append = TRUE)
    for (i in 1:dim(x)[1]) {
        cat(paste("", "node", sep = "\t"), file = file, sep = "\n", 
            append = TRUE)
        cat(paste("", "[", sep = "\t"), file = file, sep = "\n", 
            append = TRUE)
        cat(paste("\t", "id", i - 1, sep = "\t"), file = file, 
            sep = "\n", append = TRUE)
        cat(paste("\t", "label", paste("\"", dimnames(x)[[1]][i], 
            "\"", sep = ""), sep = "\t"), file = file, sep = "\n", 
            append = TRUE)
        cat(paste("\t", "graphics", sep = "\t"), file = file, 
            sep = "\n", append = TRUE)
        cat(paste("\t", "[", sep = "\t"), file = file, sep = "\n", 
            append = TRUE)
        cat(paste("\t\t", "x", stats::runif(1) * 10, sep = "\t"), file = file, 
            sep = "\n", append = TRUE)
        cat(paste("\t\t", "y", stats::runif(1) * 10, sep = "\t"), file = file, 
            sep = "\n", append = TRUE)
        cat(paste("\t\t", "type", "\"ellipse\"", sep = "\t"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("\t\t", "fill", "\"#3399FF\"", sep = "\t"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("\t\t", "outline", "\"#000000\"", sep = "\t"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("\t", "]", sep = "\t"), file = file, sep = "\n", 
            append = TRUE)
        cat(paste("\t", "LabelGraphics", sep = "\t"), file = file, 
            sep = "\n", append = TRUE)
        cat(paste("\t", "[", sep = "\t"), file = file, sep = "\n", 
            append = TRUE)
        cat(paste("\t\t", "text", paste("\"", dimnames(x)[[1]][i], 
            "\"", sep = ""), sep = "\t"), file = file, sep = "\n", 
            append = TRUE)
        cat(paste("\t", "]", sep = "\t"), file = file, sep = "\n", 
            append = TRUE)
        cat(paste("", "]", sep = "\t"), file = file, sep = "\n", 
            append = TRUE)
    }
    pat <- c("dotted", "line", "dashed")
    for (j in 1:dim(x)[3]) {
        tmp <- transf(x[, , j], lb2lb = FALSE, ord = dim(x)[1])
        for (k in 1:length(tmp)) {
            cat(paste("", "edge", sep = "\t"), file = file, sep = "\n", 
                append = TRUE)
            cat(paste("", "[", sep = "\t"), file = file, sep = "\n", 
                append = TRUE)
            cat(paste("\t", "source", as.numeric(strsplit(tmp[k], 
                ", ")[[1]])[1] - 1, sep = "\t"), file = file, 
                sep = "\n", append = TRUE)
            cat(paste("\t", "target", as.numeric(strsplit(tmp[k], 
                ", ")[[1]])[2] - 1, sep = "\t"), file = file, 
                sep = "\n", append = TRUE)
            cat(paste("\t", "graphics", sep = "\t"), file = file, 
                sep = "\n", append = TRUE)
            cat(paste("\t", "[", sep = "\t"), file = file, sep = "\n", 
                append = TRUE)
            cat(paste("\t\t", "width", 1, sep = "\t"), file = file, 
                sep = "\n", append = TRUE)
            cat(paste("\t\t", "style", paste("\"", pat[(j%%3) + 
                1], "\"", sep = ""), sep = "\t"), file = file, 
                sep = "\n", append = TRUE)
            cat(paste("\t\t", "fill", "\"#000000\"", sep = "\t"), 
                file = file, sep = "\n", append = TRUE)
            cat(paste("\t\t", "targetArrow", "\"standard\"", 
                sep = "\t"), file = file, sep = "\n", append = TRUE)
            cat(paste("\t", "]", sep = "\t"), file = file, sep = "\n", 
                append = TRUE)
            cat(paste("", "]", sep = "\t"), file = file, sep = "\n", 
                append = TRUE)
        }
    }
    cat("]", file = file, sep = "\n", append = TRUE)
}
