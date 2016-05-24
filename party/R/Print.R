
# $Id: Print.R 282 2006-10-09 14:23:30Z hothorn $

prettysplit <- function(x, inames = NULL, ilevels = NULL) {
    if (length(x) == 4)
        names(x) <- c("variableID", "ordered", "splitpoint", "splitstatistics")
    if (length(x) == 5)
        names(x) <- c("variableID", "ordered", "splitpoint", "splitstatistics",
                      "toleft")
    if (length(x) == 6)
        names(x) <- c("variableID", "ordered", "splitpoint", "splitstatistics",
                      "toleft", "table")
    if (x$ordered) {
        class(x) <- "orderedSplit"
    } else {
        class(x) <- "nominalSplit"
    }
    if (!is.null(ilevels)) {
        if (!is.null(ilevels[x[["variableID"]]]))
            attr(x$splitpoint, "levels") <- ilevels[[x[["variableID"]]]]
    }
    if (!is.null(inames)) x$variableName <- inames[x[["variableID"]]]
    return(x)
}

prettytree <- function(x, inames = NULL, ilevels = NULL) {
    names(x) <- c("nodeID", "weights", "criterion", "terminal",
                  "psplit", "ssplits", "prediction", "left", "right")
    if (is.null(inames) && extends(class(x), "BinaryTree"))
        inames <- names(x@data@get("input"))
    names(x$criterion) <- c("statistic", "criterion", "maxcriterion")
    names(x$criterion$criterion) <- inames
    names(x$criterion$statistic) <- inames

    if (x$terminal) {
        class(x) <- "TerminalNode"
        return(x)
    }

    x$psplit <- prettysplit(x$psplit, inames = inames, ilevels = ilevels)
    if (length(x$ssplit) > 0)
        x$ssplit <- lapply(x$ssplit, prettysplit, inames = inames, 
                           ilevels = ilevels)

    class(x) <- "SplittingNode"
    x$left <- prettytree(x$left, inames = inames, ilevels = ilevels)   
    x$right <- prettytree(x$right, inames = inames, ilevels = ilevels)    
    return(x)
}
 
print.TerminalNode <- function(x, n = 1, ...) {
    cat(paste(paste(rep(" ", n - 1), collapse = ""), x$nodeID, ")* ", 
                    sep = "", collapse = ""),
        "weights =", sum(x$weights), "\n")
}
 
print.SplittingNode <- function(x, n = 1, ...) {
    cat(paste(paste(rep(" ", n - 1), collapse = ""), x$nodeID, ") ", sep=""))
    print(x$psplit, left = TRUE)
    cat(paste("; criterion = ", round(x$criterion$maxcriterion, 3), 
              ", statistic = ", round(max(x$criterion$statistic), 3), "\n", 
              collapse = "", sep = ""))
    print(x$left, n + 2)
    cat(paste(paste(rep(" ", n - 1), collapse = ""), x$nodeID, ") ", sep=""))
    print(x$psplit, left = FALSE)
    cat("\n")
    print(x$right, n + 2)
}

print.orderedSplit <- function(x, left = TRUE, ...) {
    if (!is.null(attr(x$splitpoint, "levels"))) {
        sp <- attr(x$splitpoint, "levels")[x$splitpoint]
    } else {
        sp <- x$splitpoint
    }
    if (!is.null(x$toleft)) left <- as.logical(x$toleft) == left
    if (left) {
        cat(x$variableName, "<=", sp)
    } else {
        cat(x$variableName, ">", sp)
    }
}

print.nominalSplit <- function(x, left = TRUE, ...) {

    levels <- attr(x$splitpoint, "levels")

    ### is > 0 for levels available in this node
    tab <- x$table

    if (left) {
        lev <- levels[as.logical(x$splitpoint) & (tab > 0)]
    } else {
        lev <- levels[!as.logical(x$splitpoint) & (tab > 0)]
    }

    txt <- paste("{", paste(lev, collapse = ", "), "}", collapse = "", sep = "")
    cat(x$variableName, "==", txt)
}


print.BinaryTreePartition <- function(x, ...)
    print(x@tree)

print.BinaryTree <- function(x, ...) {
    cat("\n")
    cat("\t Conditional inference tree with", length(unique(where(x))), 
        "terminal nodes\n\n")
    y <- x@responses
    if (y@ninputs > 1) {
        cat("Responses:", paste(names(y@variables), 
                                collapse = ", "), "\n")
    }  else {
        cat("Response: ", names(y@variables), "\n")
    }
    inames <- names(x@data@get("input"))
    if (length(inames) > 1) {
        cat("Inputs: ", paste(inames, collapse = ", "), "\n")
    } else {
        cat("Input: ", inames, "\n")
    }
    cat("Number of observations: ", x@responses@nobs, "\n\n")
    print(x@tree)
}

print.RandomForest <- function(x, ...) {
    cat("\n")
    cat("\t Random Forest using Conditional Inference Trees\n")
    cat("\n")
    cat("Number of trees: ", length(x@ensemble), "\n")
    cat("\n")
    y <- x@responses
    if (y@ninputs > 1) {
        cat("Responses:", paste(names(y@variables),
                                collapse = ", "), "\n")
    }  else {
        cat("Response: ", names(y@variables), "\n")
    }
    inames <- names(x@data@get("input"))
    if (length(inames) > 1) {
        cat("Inputs: ", paste(inames, collapse = ", "), "\n")
    } else {
        cat("Input: ", inames, "\n")
    }
    cat("Number of observations: ", x@responses@nobs, "\n\n")
    invisible(x)
}

setMethod("show", "BinaryTree", function(object) print(object))
setMethod("show", "RandomForest", function(object) print(object))

