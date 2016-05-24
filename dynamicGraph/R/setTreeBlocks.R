"setTreeBlocks" <-
function (block.tree, vertices, root.label = "Root", N = 3, delta = ifelse(overlaying, 
    1, 0), Delta = ifelse(overlaying, 0, 1.5), d = 5, f = 1/4, 
    blockColors = NULL, overlaying = FALSE) 
{
    "positionVerticesInTreeBlocks" <- function(block.vertices, 
        vertices, stratum, A = 0, B = 100, X = 0, Y = 100, u = 0, 
        v = 100, N = 3, horizontal = FALSE) {
        "setVertex" <- function(j, m, set.positions = TRUE) {
            if (set.positions) {
                if (m < 3) 
                  if (horizontal) 
                    position <- c(A + w * (j - 0.5), (X + Y)/2)
                  else position <- c((A + B)/2, X + w * (j - 
                    0.5))
                else position <- c((A + B)/2 + 0.9 * (B - A)/2 * 
                  cos(2 * pi * j/m), (X + Y)/2 + 0.9 * (Y - X)/2 * 
                  sin(2 * pi * j/m))
                position <- c(position, rep(50, N - 2)) - rep(50, 
                  N)
                position(vertices[[block.vertices[j]]]) <<- position
            }
            stratum(vertices[[block.vertices[j]]]) <<- stratum
            blockindex(vertices[[block.vertices[j]]]) <<- stratum
        }
        m <- length(block.vertices)
        if (!is.numeric(block.vertices)) 
            block.vertices <- match(block.vertices, Names(vertices))
        if (horizontal) 
            w <- (B - A)/m
        else w <- (Y - X)/m
        for (j in 1:m) setVertex(j, m)
        return(vertices)
    }
    "blockTreeDescendants" <- function(tree) {
        result <- tree$block@stratum
        if (!is.null((tree$sub.blocks))) 
            for (i in 1:length(tree$sub.blocks)) result <- c(result, 
                blockTreeDescendants(tree$sub.blocks[[i]]))
        tree$block@descendants <- result
        return(result)
    }
    "blockChildren" <- function(sub.blocks) {
        result <- NULL
        if (!is.null((sub.blocks))) 
            for (i in 1:length(sub.blocks)) result <- c(result, 
                sub.blocks[[i]]$block@stratum)
        return(result)
    }
    "returnTreeBlocks" <- function(block.tree, label, ancestors = NULL, 
        a = 2, b = 98, x = 2, y = 98, u = 2, v = 98, N = 3, delta = ifelse(overlaying, 
            1, 0), Delta = ifelse(overlaying, 0, 2), d = 5, f = 1/3, 
        overlaying = FALSE) {
        labels <- names(block.tree)
        closed <- FALSE
        if (!is.null(block.tree$closed)) 
            closed <- block.tree$closed
        horizontal <- TRUE
        if (!is.null(block.tree$horizontal)) 
            horizontal <- block.tree$horizontal
        if (!is.null(block.tree$d) && is.numeric(block.tree$d) && 
            (length(block.tree$d) == 1)) 
            d <- block.tree$d
        if (!is.null(block.tree$f) && is.numeric(block.tree$f) && 
            (length(block.tree$f) == 1)) 
            f <- block.tree$f
        g <- NULL
        if (!is.null(block.tree$g) && is.numeric(block.tree$g) && 
            (length(block.tree$g) == 1)) 
            g <- block.tree$g
        F <- 1
        G <- NULL
        if (!is.null(block.tree$common.children)) {
            F <- 1/2
            if (!is.null(block.tree$F) && is.numeric(block.tree$F) && 
                (length(block.tree$F) == 1)) 
                F <- block.tree$F
            if (!is.null(block.tree$G) && is.numeric(block.tree$G) && 
                (length(block.tree$G) == 1)) 
                G <- block.tree$G
        }
        if (!is.null(block.tree$label) && is.character(block.tree$label) && 
            (length(block.tree$label) == 1)) 
            label <- block.tree$label
        vertices.last <- FALSE
        if (!is.null(block.tree$vertices.last) && is.logical(block.tree$vertices.last) && 
            (length(block.tree$vertices.last) == 1)) 
            vertices.last <- block.tree$vertices.last
        if (!overlaying) 
            vertices.last <- FALSE
        n <- length(block.tree)
        count.blocks <- 0
        Label <- ""
        for (i in 1:n) if (is.list(block.tree[[i]]) && (names(block.tree)[i] != 
            "common.children")) {
            count.blocks <- count.blocks + 1
            if (label == "") 
                Label <- paste(Label, names(block.tree)[i], sep = ifelse(Label == 
                  "", "", "-"))
        }
        if (label == "") 
            label <- Label
        k <- x + d + 2 * ifelse(overlaying, delta, Delta)
        if (is.null(g)) 
            if (horizontal) 
                g <- f * (y - k)
            else g <- f * (b - a)
        Stratum <<- Stratum + 1
        next.ancestors <- c(ancestors[ancestors != 0], Stratum)
        s <- 1
        i <- which(names(block.tree) == "Vertices")
        if (length(i) == 0) 
            i <- 1
        count.vertices <- 0
        A <- a
        B <- b
        X <- k
        Y <- y
        if (count.blocks > 0) 
            if (vertices.last) 
                if (horizontal) 
                  X <- y - g + Delta
                else A <- B - g + Delta
            else if (horizontal) 
                Y <- k + g - Delta
            else B <- a + g - Delta
        if ((is.numeric(block.tree[[i]]) || !is.list(block.tree[[i]])) && 
            (names(block.tree)[i] != "d") && (names(block.tree)[i] != 
            "f") && (names(block.tree)[i] != "F") && (names(block.tree)[i] != 
            "g")) {
            count.vertices <- length(block.tree[[i]])
            if (count.vertices > 0) {
                block.vertices <- block.tree[[i]]
                if (!is.numeric(block.vertices)) 
                  block.vertices <- match(block.vertices, Names(vertices))
                vertices <<- positionVerticesInTreeBlocks(block.vertices, 
                  vertices, Stratum, A, B, X, Y, u, v, horizontal = horizontal)
                if (label == "") 
                  label <- paste(Names(vertices)[block.vertices], 
                    sep = "+")
            }
            s <- 2
        }
        g <- ifelse(count.vertices == 0, 0, g)
        if (vertices.last) {
            g.f <- 0
            g.l <- g
        }
        else {
            g.f <- g
            g.l <- 0
        }
        if (is.null(G)) 
            if (horizontal) 
                G <- F * (y - k - g)
            else G <- F * (b - a - g)
        A <- a
        X <- x
        if ((overlaying || (count.blocks == 0))) {
            B <- b
            Y <- y
        }
        else if (horizontal) 
            Y <- k + g.f - Delta
        else B <- a + g.f - Delta
        position <- matrix(c(c(A, X, rep(u, N - 2)), c(B, Y, 
            rep(v, N - 2))), nrow = 2, byrow = T) - rep(50, N)
        depth <- length(ancestors)
        depth <- 2 * depth + 1
        if (is.element(0, ancestors)) 
            depth <- depth - 1
        block.color <- "Grey"
        if (!is.null(blockColors) && (length(blockColors) >= 
            depth) && is.character(blockColors[depth])) 
            block.color <- blockColors[depth]
        else {
            if (!is.null(blockColors) && (length(blockColors) == 
                1) && is.numeric(blockColors[1])) 
                depth <- round(depth * blockColors[1])
            block.color <- paste("grey", 100 - depth, sep = "")
        }
        Block <- new("dg.Block", stratum = Stratum, index = -Stratum, 
            position = position, color = block.color, closed = closed, 
            label = label, ancestors = if (is.null(ancestors)) 
                0
            else ancestors)
        result <- NULL
        if (count.blocks > 0) {
            if (horizontal) 
                w <- (b - a)/count.blocks
            else {
                if (overlaying) 
                  if ((x + 15 < y)) 
                    x <- k
                w <- (y - x)/count.blocks
            }
            if (horizontal) {
                X <- k + g.f
                Y <- X + G - delta
            }
            else {
                A <- a + max(delta, g.f)
                B <- A + G - delta
            }
            j <- -1
            for (i in s:n) if (is.list(block.tree[[i]]) && (names(block.tree)[i] != 
                "common.children")) {
                j <- j + 1
                if (horizontal) {
                  A <- a + delta + w * j
                  B <- a - delta + w * (1 + j)
                }
                else {
                  X <- x + delta + w * j
                  Y <- x - delta + w * (1 + j)
                }
                if ((!overlaying) && (j + 1 < count.blocks)) 
                  if (horizontal) 
                    B <- B - Delta
                  else Y <- Y - Delta
                result <- append(result, list(returnTreeBlocks(block.tree[[i]], 
                  labels[i], next.ancestors, A, B, X, Y, u, v, 
                  N = N, delta = delta, Delta = Delta, d = d, 
                  f = f, overlaying = overlaying)))
            }
            if (horizontal) {
                A <- a + delta
                B <- b - delta
                X <- X + delta + G + Delta
                Y <- y - delta - g.l
            }
            else {
                X <- x + delta
                Y <- y - delta
                A <- A + delta + G + Delta
                B <- b - delta - g.l
            }
            if (!is.null(block.tree$common.children)) 
                result <- append(result, list(returnTreeBlocks(block.tree$common.children, 
                  "", c(0, next.ancestors), A, B, X, Y, u, v, 
                  N = N, delta = delta, Delta = Delta, d = d, 
                  f = f, overlaying = overlaying)))
            des <- blockTreeDescendants(list(block = Block, sub.blocks = result))
            Block@descendants <- des[!((des == Block@stratum) | 
                (des == 0))]
            cld <- blockChildren(result)
            Block@children <- cld[!((cld == Block@stratum) | 
                (cld == 0))]
        }
        return(list(block = Block, sub.blocks = result))
    }
    for (i in seq(along = vertices)) position(vertices[[i]]) <- c(100 * 
        (i - 0.5)/length(vertices) - 50, -48.5, rep(0, N - 2))
    Stratum <- 0
    result <- returnTreeBlocks(block.tree, label = root.label, 
        NULL, N = N, delta = delta, Delta = Delta, d = d, f = f, 
        overlaying = overlaying)
    return(list(Vertices = vertices, BlockTree = result))
}
