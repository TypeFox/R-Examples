countsToBinomial <- function(xtab) {
    ## make square if necessary
    if (nrow(xtab) != ncol(xtab) || rownames(xtab) != colnames(xtab)) {
        dat <- as.data.frame(xtab)
        lev <- union(rownames(xtab), colnames(xtab))
        dat[,1] <- factor(dat[,1], levels = lev)
        dat[,2] <- factor(dat[,2], levels = lev)
        xtab <- tapply(dat[,3], dat[1:2], sum)
        xtab[is.na(xtab)] <- 0
    }
    ##assumes square
    players <- rownames(xtab)
    comb <- combinations(nrow(xtab), 2)
    won <- xtab[comb]
    lost <- t(xtab)[comb]
    res <- !(won == 0 & lost == 0)
    player1 <- factor(players[comb[,1]], levels = players)[res]
    player2 <- factor(players[comb[,2]], levels = players)[res]
    data.frame(player1, player2, win1 = won[res], win2 = lost[res])
}

