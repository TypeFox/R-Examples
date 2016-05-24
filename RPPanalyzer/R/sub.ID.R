`sub.ID` <-
function(x) {
        strings <- rownames(x[[1]])
        id <-   sub("^[A-Za-z]*?-[0-9]*?-[0-9]*?-0","",strings)
        ident <- sub("-", "", id)
        rownames (x[[1]]) <- ident
        rownames (x[[2]]) <- ident
        rownames (x[[3]]) <- ident
        return (x)
}

