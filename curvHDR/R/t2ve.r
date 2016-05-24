
t2ve <- function (triangles)
{
    vb <- rbind(triangles$v1, triangles$v2, triangles$v3)
    vbmin <- min(vb)
    vbmax <- max(vb)
    S <- 10^5
    score <- function(v, d) floor(as.vector(v %*% d))
    scale <- function(v) (1 - 1/S) * (v - vbmin)/(vbmax - vbmin)
    d <- c(S, S^2, S^3)
    scores <- score(scale(vb), d)
    vb <- vb[!duplicated(scores), ]
    scores <- score(scale(vb), d)
    ib <- rbind(match(score(scale(triangles$v1), d), scores), 
        match(score(scale(triangles$v2), d), scores), match(score(scale(triangles$v3), 
            d), scores))
    return(list(vb = t(vb), ib = ib))
}


