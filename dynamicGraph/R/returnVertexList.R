"returnVertexList" <-
function (names, labels = NULL, types = NULL, strata = NULL, 
    line = FALSE, N = 3, colors = ifelse(types == "TextVertex", 
        "FloralWhite", "DarkRed"), vertexClasses = validVertexClasses()) 
{
    "newNodeList" <- function(list) return(new("dg.NodeList", 
        nodeList = list))
    "newVertexList" <- function(list) return(new("dg.VertexList", 
        nodeList = list))
    n <- length(names)
    if (length(colors) == 1) 
        colors <- rep(colors, n)
    if (length(colors) == 0) 
        colors <- rep("DarkRed", n)
    if (!is.null(labels) && !(length(labels) == n)) 
        warning("Invalid length of argument labels")
    if (!is.null(types) && (length(types) == 1)) 
        types <- rep(types, n)
    if (!is.null(types) && !(length(types) == n)) 
        warning("Invalid length of argument types")
    if (!is.null(strata) && !(length(strata) == n)) 
        warning("Invalid length of argument strata")
    if (!is.null(colors) && !(length(colors) == n)) 
        warning("Invalid length of argument colors")
    result <- vector("list", n)
    for (i in seq(along = names)) {
        if (line) 
            position <- c(-45, (40 * (1 - (2 * ((i - 1) + 0.5)/n))))
        else position <- c(40 * cos(2 * pi * i/n), 40 * sin(2 * 
            pi * i/n))
        position <- c(position, rep(0, max(0, N - 2)))
        if (is.null(labels)) 
            label <- names[i]
        else label <- labels[i]
        if (is.null(types)) 
            type <- vertexClasses[, 1][1]
        else type <- types[i]
        if (is.null(strata)) 
            stratum <- 0
        else stratum <- strata[i]
        if (!is.na(type) && type == "TextVertex") 
            prototype <- "dg.TextVertex"
        else {
            prototype <- "dg.Vertex"
            x <- match(type, vertexClasses[, 1])
            if (!is.null(x) && !all(is.na(x))) 
                prototype <- paste(vertexClasses[, 2][x])
        }
        if (prototype == "dg.TextVertex") 
            index <- -i
        else index <- i
        result[[i]] <- new(prototype, name = names[i], label = label, 
            index = index, position = position, blockindex = 0, 
            stratum = stratum, color = colors[i])
    }
    class(result) <- "dg.VertexList"
    names(result) <- Names(result)
    return(result)
}
