`mantel.pertables` <-
function (pertab, env, dist.method = "bray", binary = FALSE, 
    cor.method = "pearson", permutations = 100) 
{
    require(vegan)
    mantel.test <- function(z, env) {
        mantel.st <- mantel(vegdist(z, method = dist.method, 
            binary = binary), dist(env), method = cor.method, 
            permutations = permutations)
        mantel.r <- -mantel.st$statistic
        mantel.p <- mantel.st$signif
        return(c(mantel.r, mantel.p))
    }
    results <- sapply(pertab$pertables, function(x) mantel.test(x, 
        env))
    row.names(results) <- c("r", "p-value")
    mantel.quant <- apply(results, 1, quantile, c(0, 0.005, 0.025, 
        0.5, 0.975, 0.995, 1))
    vdper <- lapply(pertab$pertables, function(x) 1 - vegdist(x, 
        method = dist.method, binary = binary))
    z <- pertab$raw
    mantel.raw <- mantel(vegdist(z, method = dist.method, binary = binary), 
        dist(env), method = cor.method, permutations = permutations)
    mantel.r <- -mantel.raw$statistic
    ptax <- ((rank(c(mantel.r, results[1, ])))/(length(results[1, 
        ]) + 1))[1]
    ptax <- ifelse(ptax <= 0.5, ptax, 1 - ptax)
    vd <- 1 - vegdist(pertab$raw, method = dist.method, binary = binary)
    env.dist <- as.vector(dist(env))
    mantel.output <- list(mantel = list(mantel.raw = mantel.raw, 
        ptax = ptax), simulation = list(results = results, mantel.quant = mantel.quant, 
        vegdist = vdper), raw = list(vegdist = vd, env.dist = env.dist))
    class(mantel.output) <- c("mantel.pertables", class(mantel.output))
    return(mantel.output)
}

