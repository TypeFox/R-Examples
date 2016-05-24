Fbwidths.by.x <-
function (tab.x, tab.xy, tab.xz, compress.sum=FALSE) 
{
    N <- sum(tab.xy) + sum(tab.xz)
    prop.x <- prop.table(tab.x)
    prop.xy <- prop.table(tab.xy)
    prop.xz <- prop.table(tab.xz)
    
    lab.x <- names(dimnames(tab.x))
    if (all(nchar(lab.x) == 0)) 
        lab.x <- paste("x", 1:length(lab.x), sep = "")
    names(attr(tab.x, "dimnames")) <- lab.x
    
    lab.xy <- names(dimnames(tab.xy))
    if (all(nchar(lab.xy) == 0)) 
        lab.xy <- c(lab.x, "y")
    names(attr(tab.xy, "dimnames")) <- lab.xy
    
    lab.y <- setdiff(lab.xy, lab.x)
    p.y <- match(lab.y, lab.xy)
    
    lab.xz <- names(dimnames(tab.xz))
    if (all(nchar(lab.xz) == 0)) 
        lab.xz <- c(lab.x, "z")
    names(attr(tab.xz, "dimnames")) <- lab.xz
    
    lab.z <- setdiff(lab.xz, lab.x)
    p.z <- match(lab.z, lab.xz)
#    
    n.x <- length(lab.x)
    appo.var <- as.list(lab.x)
    for (k in 2:n.x) {
        b <- combn(lab.x, k)
        b <- data.frame(b, stringsAsFactors = FALSE)
        appo.var <- c(appo.var, as.list(b))
    }
    
    H <- length(appo.var)
    out.rng <- as.list(as.numeric(H))
    av.rng <- matrix(NA, H, 8)
#    av.rng <- matrix(NA, H, 9)
#    all.H <- matrix(NA, H, 5)
#    all.U <- matrix(NA, H, 2)
    
    for (h in 1:H) {
        lab <- appo.var[[h]]
        
        p.x <- match(lab, lab.x)
       
        xx <- margin.table(prop.x, p.x)
        av.rng[h, 1] <- length(xx)
        av.rng[h, 2] <- sum(xx == 0)
        p.xy <- match(c(lab, lab.y), lab.xy)
        xy <- margin.table(prop.xy, p.xy)
        av.rng[h, 3] <- length(xy)
        av.rng[h, 4] <- sum(xy == 0)
        p.xz <- match(c(lab, lab.z), lab.xz)
        xz <- margin.table(prop.xz, p.xz)
        av.rng[h, 5] <- length(xz)
        av.rng[h, 6] <- sum(xz == 0)

        fb <- Frechet.bounds.cat(xx, xy, xz, print.f = "tables")
        appo <- data.frame(fb$low.cx)
        out.rng[[h]] <- data.frame(appo[, 1:2], lower = c(fb$low.cx), 
                                   upper = c(fb$up.cx), width = c(fb$up.cx - fb$low.cx))
        
        av.rng[h, 7] <- fb$uncertainty[2]
        av.rng[h, 8] <- fb$uncertainty[2] / fb$uncertainty[1]
#        av.rng[h, 9] <- fb$uncertainty[3]
#        all.H[h, ] <- fb$H
#        all.U[h, ] <- fb$U
    }
    lab.list <- paste("|", lapply(appo.var, paste, collapse = "+"), 
                      sep = "")
    n.vars <- lapply(appo.var, length)
    av.rng <- data.frame(x.vars = unlist(n.vars), 
                         x.cells = av.rng[, 1], x.freq0 = av.rng[, 2], 
                         xy.cells = av.rng[, 3], xy.freq0 = av.rng[, 4],
                         xz.cells = av.rng[, 5], xz.freq0 = av.rng[, 6],
                         av.width = av.rng[, 7], rel.av.width = av.rng[, 8])
#                         delta.CMS=av.rng[, 9])
    row.names(av.rng) <- paste("|", lapply(appo.var, paste, collapse = "*"), 
                               sep = "")
    av.rng.0 <- c(x.vars=0, x.cells=NA, x.freq0=NA, 
                  xy.cells = NA, xy.freq0 = NA,
                  xz.cells = NA, xz.freq0 = NA,
                  av.width = fb$uncertainty[1], rel.av.width = 1)
    av.rng <- rbind(unconditioned=av.rng.0, av.rng)
#    colnames(all.H) <- names(fb$H)
#    colnames(all.U) <- names(fb$U)
#    row.names(all.H) <- rownames(all.U) <- paste("|", lapply(appo.var, paste, collapse = "+"), sep = "")
    aa <- n.x - av.rng$x.vars
    ord.lab <- order(aa, av.rng$av.width, decreasing = TRUE)
    av.rng <- av.rng[ord.lab, ]
    if(compress.sum){
        sp.av <- split(av.rng, av.rng$x.vars)
        G <- length(sp.av)
        sp.new <- as.list(G)
        sp.new[[1]] <- sp.av[[1]]
        sp.new[[2]] <- sp.av[[2]]
        for(g in 3:G){
            min.p <- min(sp.av[[(g-1)]][,"av.width"])
            tst <- sp.av[[g]][,"av.width"] <= min.p
            sp.new[[g]] <- sp.av[[g]][tst,]
        }
        av.rng <- do.call("rbind", sp.new)
    }
    out.rng[[(H + 1)]] <- av.rng
    names(out.rng) <- c(lab.list, "sum.unc")
#    out.rng$all.H <- all.H[ord.lab,]
#    out.rng$all.U <- all.U[ord.lab,]

    out.rng
}