# Automatically generated from all.nw using noweb
print.pedigree <- function(x, ...) {
    cat("Pedigree object with", length(x$id), "subjects")
    if (!is.null(x$famid)) cat(", family id=", x$famid[1], "\n")
    else cat("\n")
    cat("Bit size=", bitSize(x)$bitSize, "\n")
    }

print.pedigreeList <- function(x, ...) {
    cat("Pedigree list with", length(x$id), "total subjects in",
        length(unique(x$famid)), "families\n")
    }
