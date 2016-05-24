"setBlocks" <-
function (block.list, vertices, labels = NULL, right.to.left = FALSE, 
    nested.blocks = FALSE, blockColors = NULL, color = "Grey", N = 3) 
{
    "newBlockList" <- function(list) return(new("dg.BlockList", 
        nodeList = list))
    "returnDefaultBlocks" <- function(strata = NULL, labels = NULL, 
        n = length(strata), right.to.left = FALSE, nested.blocks = FALSE, 
        visible = TRUE, N = 3) {
        result <- vector("list", n)
        A <- 2/3
        B <- 3/4
        C <- 1/2
        D <- c(A, B, rep(C, N - 2))
        for (j in 1:n) {
            i <- strata[j]
            if (nested.blocks) {
                position <- matrix(95 * (c(D * (i - 1)/n,
                                   1 - (1 - D) * (i - 1)/n) - 0.5), 
                                   nrow = 2, byrow = TRUE)
            }
            else position <- matrix(c(c(
                80/n * (i - 1) + 20/(1 + n) * i - 50, -45, rep(-50, N - 2)), 
                c((80/n + 20/(1 + n)) * i - 50, 45, rep(50, N - 2))), 
                nrow = 2, byrow = TRUE)
            if (is.null(labels)) 
                label <- paste("Block", j)
            else label <- labels[i]
            if (!is.null(blockColors) && (length(blockColors) >= j)) 
                color <- blockColors[i]
            result[[i]] <- new("dg.Block", stratum = j, index = -j, 
                position = position, closed = FALSE, visible = TRUE, 
                color = color, label = label, parent = 0, ancestors = 0)
        }
        names(result) <- labels
        return(result)
    }
    "positionVerticesInDefaultBlocks" <- function(block.list, 
        vertices, right.to.left = FALSE, nested.blocks = FALSE, 
        n = length(block.list), draw.blocks = FALSE, N = N) {
        "setVertex" <- function(i, j, n, m, right.to.left, 
                                set.positions = TRUE) {
            if (set.positions) {
                if (nested.blocks) {
                  A <- 2/3
                  B <- 3/4
                  C <- 1/2
                  D <- c(A, B, rep(C, N - 2))
                  posA <- matrix(
                    95 * (c(D * (i - 1)/n, 1 - (1 - D) * (i - 1)/n) - 0.5), 
                    nrow = 2, byrow = TRUE)
                  posB <- matrix(
                    95 * (c(D * i/n, 1 - (1 - D) * i/n) - 0.5), 
                    nrow = 2, byrow = TRUE)
                  pos <- (posA + posB)/2
                  position <- pos[1, ]
                  if (m > 1) 
                    position[2] <- position[2] + 
                      (pos[2, 2] - pos[1, 2]) * (j - 1)/(m - 1)
                  position[3] <- 0
                }
                else if (m < 3) 
                  position <- c(
                    80/n * (i - 0.5) + 20/(1 + n) * i - 50, 
                    80/m * (j - 0.5) + 10/(1 + m) * j - 45, 
                    rep(0, N - 2))
                else position <- c(
                    80/n * (i - 0.5) + 20/(1 + n) * i + 
                       0.8 * 40 * cos(2 * pi * j/m)/n - 50, 
                    0 + 40 * sin(2 * pi * j/m), 
                    rep(0 + 40 * cos(2 * pi * j/m), N - 2))
                position(vertices[[name.list[j]]]) <<- position
            }
            stratum(vertices[[name.list[j]]]) <<- ifelse(right.to.left, 
                n - i, i)
            blockindex(vertices[[name.list[j]]]) <<- stratum(
                vertices[[name.list[j]]])
        }
        for (i in seq(along = vertices)) position(vertices[[i]]) <- c(100 * 
            (i - 0.5)/length(vertices) - 50, -48.5, rep(0, N - 2))
        if (right.to.left) 
            i.list <- rev(1:n)
        else i.list <- 1:n
        for (ii in 1:n) {
            name.list <- block.list[[ii]]
            i <- i.list[ii]
            m <- length(name.list)
            if (!is.numeric(name.list)) 
                name.list <- nameToVertexIndex(name.list, vertices)
            for (j in 1:m) {
                setVertex(i, j, n, m, right.to.left, draw.blocks)
            }
        }
        return(vertices)
    }
    labels <- names(block.list)
    Blocks <- returnDefaultBlocks(strata = 1:length(block.list), 
        labels = labels, n = length(block.list), 
        right.to.left = right.to.left, 
        nested.blocks = nested.blocks, visible = TRUE, N = N)
    Vertices <- positionVerticesInDefaultBlocks(block.list, vertices, 
        right.to.left = right.to.left, nested.blocks = nested.blocks, 
        n = length(block.list), draw.blocks = TRUE, N = N)
    class(Blocks) <- "dg.BlockList"
    return(list(Blocks = Blocks, Vertices = Vertices))
}
