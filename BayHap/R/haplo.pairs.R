 haplo.pairs<-function (geno) 
{
    i <- 1
    nlocus <- ncol(geno)/2
    H <- haplo.list.all(floor(nlocus))
    haplos.all <- apply(H, 1, FUN = paste, collapse = "")
    long <- NULL
    cod.haplo.cromo.total <- list()
    na <- NULL
    for (i in 1:nrow(geno)) {
        cod.haplo.cromo1 <- NULL
        cod.haplo.cromo2 <- NULL
        if (any(geno[i, ] == 8)) {
            j <- 1
            haplo1.ini <- NULL
            haplo2.ini <- NULL
            mat <- NULL
            for (j in 1:ncol(geno)) {
                ifelse(j%%2 == 0, haplo2.ini <- cbind(haplo2.ini, 
                  geno[i, j]), haplo1.ini <- cbind(haplo1.ini, 
                  geno[i, j]))
            }
            pos.na <- seq(1:nlocus)[haplo1.ini == 8]
            n.na <- sum(haplo1.ini == 8)
            h.na <- haplo.list.all(n.na)
            k <- 1
            l <- 1
            haplo1.expand <- NULL
            haplo2.expand <- NULL
            for (l in 1:nrow(h.na)) {
                haplo1.expand <- rbind(haplo1.expand, haplo1.ini)
                haplo2.expand <- rbind(haplo2.expand, haplo2.ini)
            }
            haplo1.expand[, pos.na] <- h.na
            haplo2.expand[, pos.na] <- h.na
            a <- 1
            row.haplo2 <- 0
            for (a in 1:nrow(haplo1.expand)) {
                haplo1 <- haplo1.expand[a, ]
                row.haplo2 <- row.haplo2 + 1
                b <- row.haplo2
                for (b in row.haplo2:nrow(haplo2.expand)) {
                  haplo2 <- haplo2.expand[b, ]
                  hetero <- rep(0, nlocus)
                  hetero[haplo1 != haplo2] <- 1
                  nhetero <- sum(hetero)
                  ifelse(nhetero == 0, rows <- 1, rows <- 2^(nhetero - 
                    1))
                  cromos <- haplo.list(nlocus, haplo1, haplo2, 
                    hetero, nhetero, rows, nlocus)
                  l <- 1
                  m <- 1
                  haplos.indiv <- matrix(rep(9, 2 * nlocus * 
                    rows), nrow = 2 * rows)
                  for (l in 1:rows) {
                    for (m in 1:nlocus) {
                      haplos.indiv[2 * l - 1, m] <- cromos[[1]][l + 
                        rows * (m - 1)]
                      haplos.indiv[2 * l, m] <- cromos[[2]][l + 
                        rows * (m - 1)]
                    }
                  }
                  haplos.indiv <- apply(haplos.indiv, 1, FUN = paste, 
                    collapse = "")
                  k <- 1
                  pos <- NULL
                  for (k in 1:(2 * rows)) {
                    pos <- cbind(pos, seq(0, (2^nlocus) - 1)[haplos.indiv[k] == 
                      haplos.all])
                  }
                  k <- 1
                  pos1 <- NULL
                  pos2 <- NULL
                  for (k in 1:(2 * rows)) {
                    if (k%%2 != 0) 
                      pos1 <- cbind(pos1, pos[k])
                    if (k%%2 == 0) 
                      pos2 <- cbind(pos2, pos[k])
                  }
                  cod.haplo.cromo1 <- c(cod.haplo.cromo1,pos1)
                  cod.haplo.cromo2 <- c(cod.haplo.cromo2,pos2)
                }
            }
            mat <- rbind(cod.haplo.cromo1, cod.haplo.cromo2)
            mat.ord <- apply(as.data.frame(mat), 2, FUN = function(x) {
                x <- x[order(x)]
            })
            cod.haplo.cromo <- apply(as.data.frame(mat.ord), 
                2, FUN = paste, collapse = "-")
            cod.haplo.cromo <- unique(cod.haplo.cromo)
            long <- c(long, length(cod.haplo.cromo))
            cod.haplo.cromo.total[[length(cod.haplo.cromo.total) + 
                1]] <- cod.haplo.cromo
        }
        else {
            j <- 1
            haplo1 <- NULL
            haplo2 <- NULL
            for (j in 1:ncol(geno)) {
                ifelse(j%%2 == 0, haplo2 <- cbind(haplo2, geno[i, 
                  j]), haplo1 <- cbind(haplo1, geno[i, j]))
            }
            hetero <- rep(0, nlocus)
            hetero[haplo1 != haplo2] <- 1
            nhetero <- sum(hetero)
            ifelse(nhetero == 0, rows <- 1, rows <- 2^(nhetero - 
                1))
            cromos <- haplo.list(nlocus, haplo1, haplo2, hetero, 
                nhetero, rows, nlocus)
            l <- 1
            m <- 1
            haplos.indiv <- matrix(rep(9, 2 * nlocus * rows), 
                nrow = 2 * rows)
            for (l in 1:rows) {
                for (m in 1:nlocus) {
                  haplos.indiv[2 * l - 1, m] <- cromos[[1]][l + 
                    rows * (m - 1)]
                  haplos.indiv[2 * l, m] <- cromos[[2]][l + rows * 
                    (m - 1)]
                }
            }
            haplos.indiv <- apply(haplos.indiv, 1, FUN = paste, 
                collapse = "")
            k <- 1
            pos <- NULL
            for (k in 1:(2 * rows)) {
                pos <- cbind(pos, seq(0, (2^nlocus) - 1)[haplos.indiv[k] == 
                  haplos.all])
            }
            k <- 1
            pos1 <- NULL
            pos2 <- NULL
            for (k in 1:(2 * rows)) {
                if (k%%2 != 0) 
                  pos1 <- cbind(pos1, pos[k])
                if (k%%2 == 0) 
                  pos2 <- cbind(pos2, pos[k])
            }
            cod.haplo.cromo1 <- c(cod.haplo.cromo1, pos1)
            cod.haplo.cromo2 <- c(cod.haplo.cromo2, pos2)
            mat <- rbind(cod.haplo.cromo1, cod.haplo.cromo2)
            mat.ord <- apply(as.data.frame(mat), 2, FUN = function(x) {
                x <- x[order(x)]
            })
            cod.haplo.cromo <- apply(as.data.frame(mat.ord), 
                2, FUN = paste, collapse = "-")
            cod.haplo.cromo <- unique(cod.haplo.cromo)
            long <- c(long, length(cod.haplo.cromo))
            cod.haplo.cromo.total[[length(cod.haplo.cromo.total) + 
                1]] <- cod.haplo.cromo
        }
    }
   return(list(long=long,cod.haplo.cromo.total=cod.haplo.cromo.total))
}
