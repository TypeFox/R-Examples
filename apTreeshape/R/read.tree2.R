read.tree2 <- function (file = "", format = "Newick", rooted = TRUE, text = NULL, 
    tree.names = NULL, skip = 0, comment.char = "#", ...) 
{
    if (!is.null(text)) {
        if (!is.character(text)) 
            stop("argument `text' must be of mode character")
        tree <- text
    }
    else {
        tree <- scan(file = file, what = character(), sep = "\n", 
            quiet = TRUE, skip = skip, comment.char = comment.char,...)
    }
    if (identical(tree, character(0))) {
        warning("empty character string.")
        return(NULL)
    }
    tree <- gsub("[ \t]", "", tree)
    tsp <- unlist(strsplit(tree, NULL))
    ind <- which(tsp == ";")
    nb.tree <- length(ind)
    x <- c(1, ind[-nb.tree] + 1)
    y <- ind - 1
    STRING <- list()

if (is.na(y[1])) {cat("Invalid tree\n"); return(NULL)}
else {
    for (i in 1:nb.tree) {

STRING[[i]] <- paste(tsp[x[i]:y[i]], 
        sep = "", collapse = "")}
    list.obj <- list()
    for (i in 1:nb.tree) {
        list.obj[[i]] <- if (length(grep(":", STRING[[i]]))) 
            tree.build(STRING[[i]])
        else clado.build(STRING[[i]])
        if (sum(list.obj[[i]]$edge[, 1] == "-1") == 1) {
            warning("The root edge is apparently not correctly represented\nin your tree: this may be due to an extra pair of\nparentheses in your file; the returned object has been\ncorrected but your file may not be in a valid Newick\nformat")
            ind <- which(list.obj[[i]]$edge[, 1] == "-1")
            list.obj[[i]]$root.edge <- list.obj[[i]]$edge.length[ind]
            list.obj[[i]]$edge.length <- list.obj[[i]]$edge.length[-ind]
            list.obj[[i]]$edge <- list.obj[[i]]$edge[-ind, ]
            for (j in 1:length(list.obj[[i]]$edge)) if (as.numeric(list.obj[[i]]$edge[j]) < 
                0) 
                list.obj[[i]]$edge[j] <- as.character(as.numeric(list.obj[[i]]$edge[j]) + 
                  1)
            if (sum(list.obj[[i]]$edge[, 1] == "-1") == 1) 
                stop("There is apparently two root edges in your file:\ncannot read tree file")
        }
    }
    if (nb.tree == 1) 
        list.obj <- list.obj[[1]]
    else {
        if (is.null(tree.names)) 
            names(list.obj) <- paste("tree", 1:nb.tree, sep = "")
        else names(list.obj) <- tree.names
        class(list.obj) <- c("multi.tree", "phylo")
    }
    list.obj
}
}