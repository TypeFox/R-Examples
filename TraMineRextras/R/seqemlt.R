
seqemlt <- function(seqdata, a = 1, b = 1, weighted = TRUE) {


    ## Définition des fonctions

    ## La fonction freq définit la fréquence de chaque situation (l'état indicé par
    ## le temps) Elle permet de suivre l'évolution de la fréquence de chaque état
    ## dans le temps
    freq <- function(seq) {

        situation.time <- rep(1:dim(seq)[2], each = length(alphabet(seq)))
        situation.states.time <- paste(alphabet(seq), situation.time, sep = ".")
        situation.states <- rep(alphabet(seq), times = dim(seq)[2])
        sit2 <- situation.time
        for (i in 1:length(sit2)) {
            sit2[i] = length(which(seq[, situation.time[i]] == situation.states[i]))
        }
        names(sit2) <- situation.states.time
        return(sit2)
    }

    ## La fonction disjonctif code les séquences sous la forme du disjonctif
    ## complet

    disjonctif <- function(seq) {

        situation.time <- rep(1:dim(seq)[2], each = length(alphabet(seq)))
        situation.states.time <- paste(alphabet(seq), situation.time, sep = ".")
        situation.states <- rep(alphabet(seq), times = dim(seq)[2])
        x <- array(0, dim = c(dim(seq)[1], length(situation.states.time)))
        for (i in 1:dim(seq)[1]) {
            for (j in 1:dim(seq)[2]) {
                x[i, which(situation.time == j & situation.states == seq[i, j])] <- 1
            }
        }
        colnames(x) <- situation.states.time
        return(x)
    }


    ## La fonction profil fournit pour chaque situation (état indicé par le temps)
    ## le vecteur profil des transitions d'une situation S vers son futur Soit la
    ## probabilite de transition de S vers touts les situations futures S' Pour tous les
    ## pas Les pas (le temps séparant S et S') étant pondéré par le fonction de
    ## décroissance du temps de 1/(param_a *t param_b) Le profil est NA si la
    ## situation n'apparait pas dans la base des séquence Elle attend en entrée la
    ## séquence sous forme alphabetique (seq) et de disjonctif complet (disj)

    transrate <- function(seq, disj, poids = NULL) {
        state <- alphabet(seq)
        situation.time <- rep(1:dim(seq)[2], each = length(alphabet(seq)))
        situation.states.time <- paste(alphabet(seq), situation.time, sep = ".")
        situation.states <- rep(alphabet(seq), times = dim(seq)[2])

        x <- array(0, dim = c(length(situation.states.time), length(situation.states.time)))
        disj.pond <- disj
        if (!is.null(poids)) {
            disj.pond <- diag(poids) %*% disj
        }

        for (i in 1:length(situation.states.time)) {
            nb <- sum(disj.pond[, i])

            for (j in i:length(situation.states.time)) {
                t <- situation.time[j] - situation.time[i]
                u = crossprod(disj.pond[, i], disj.pond[, j])
                if (nb > 0) {
                  x[i, j] <- u/nb
                } else x[i, j] <- NA

            }

            alpha <- sum(x[i, ])

            x[i, ] <- x[i, ]
        }
        colnames(x) <- situation.states.time
        rownames(x) <- situation.states.time
        return(x)
    }




    profil <- function(seq, transrate, param_a = 1, param_b = 1) {
        state <- alphabet(seq)
        situation.time <- rep(1:dim(seq)[2], each = length(alphabet(seq)))
        situation.states.time <- paste(alphabet(seq), situation.time, sep = ".")
        situation.states <- rep(alphabet(seq), times = dim(seq)[2])

        x <- array(0, dim = c(length(situation.states.time), length(situation.states.time)))

        for (i in 1:length(situation.states.time)) {

            for (j in i:length(situation.states.time)) {
                t <- situation.time[j] - situation.time[i]
                beta <- param_a * t + param_b
                if (!is.na(transrate[i, j])) {
                  x[i, j] <- (transrate[i, j]/beta)
                } else x[i, j] <- NA

            }

            alpha <- sum(x[i, ])

            x[i, ] <- x[i, ]/alpha
        }
        colnames(x) <- situation.states.time
        rownames(x) <- situation.states.time
        return(x)
    }

    ## La fonction distsquare est une fonction intermédaires.  Elle définit la
    ## distance euclidienne entre les profil La distance est NA si le profil n'est
    ## pas définit c-à-d si la situation n'apparait pas dans la base des séquence
    ## Elle attend en entrée la séquence et le profil


    distsquare <- function(seq, profil) {
        situation.time <- rep(1:dim(seq)[2], each = length(alphabet(seq)))
        situation.states.time <- paste(alphabet(seq), situation.time, sep = ".")
        situation.states.time.freq <- freq(seq)
        freq <- situation.states.time.freq[which(situation.states.time.freq != 0)]
        prof <- profil[which(situation.states.time.freq != 0), which(situation.states.time.freq !=
            0)]
        D1 <- matrix(0, ncol = dim(prof)[1], nrow = dim(prof)[1])
        profilLigneSum <- as.matrix(apply(prof, 2, sum))
        ### ici
        for (i in 1:dim(D1)[1]) {
            for (j in 1:i) {
                u <- 0
                dp <- (prof[i, ] - prof[j, ])
                u <- crossprod((dp/profilLigneSum), dp)
                D1[i, j] <- u
                D1[j, i] <- u
            }
        }
        D <- matrix(NA, ncol = dim(profil)[1], nrow = dim(profil)[1])
        colnames(D) <- situation.states.time
        rownames(D) <- situation.states.time
        D[which(situation.states.time.freq != 0), which(situation.states.time.freq !=
            0)] <- D1



        return(D)
    }




    ## La fonction benz() est une fonction intermédaires.  Elle définit la matrice
    ## des covariances entre les profil à partir de la matrice des distances Benz
    ## se limite aux situation réelle (celle qui se sont produites au moins une
    ## fois) et sert d'entrée à la fonction recode()et cor()

    benz <- function(seq, d) {
        situation.time <- rep(1:dim(seq)[2], each = length(alphabet(seq)))
        situation.states.time <- paste(alphabet(seq), situation.time, sep = "")
        situation.states.time.freq <- freq(seq)
        D <- d[which(situation.states.time.freq != 0), which(situation.states.time.freq !=
            0)]

        DiMoy <- as.matrix(apply(D, 1, mean))
        DMoyj <- t(as.matrix(apply(D, 2, mean)))
        moy <- mean(DiMoy)
        Dl <- matrix(DMoyj, ncol = dim(D)[1], nrow = dim(D)[2], byrow = TRUE)
        Dc <- matrix(DiMoy, ncol = dim(D)[1], nrow = dim(D)[2])
        Dlc <- matrix(moy, ncol = dim(D)[1], nrow = dim(D)[2])
        benz <- (-1/2) * (D - Dl - Dc + Dlc)
        return(benz)
    }

    ## La fonction cor() fournit la corrélation entre les situations (indicatrices
    ## de connaitre la situation s=e.t, soit d'être dans l'état e à l'instant t) La
    ## corrélation entre deux situations est la corrélation entre les situations
    ## dans l'espace des profils de futurs Deux situations sont fortement corrélées
    ## si leurs futurs sont proches deux situations ont des futurs proches si les
    ## taux de transition de chacunes vers les autres états et eux mêmes sont
    ## identiques et ce au cours du temps Des futurs proches à court termes génèrent
    ## des corrélations plus fortes que des futurs proches à long term (cf
    ## paramètres a et b)

    corel <- function(seq, benz) {
        situation.time <- rep(1:dim(seq)[2], each = length(alphabet(seq)))
        situation.states.time <- paste(alphabet(seq), situation.time, sep = "")
        situation.states.time.freq <- freq(seq)

        x <- cor(benz, method = "pearson")


        cor <- matrix(NA, ncol = length(situation.states.time), nrow = length(situation.states.time))
        cor[which(situation.states.time.freq != 0), which(situation.states.time.freq !=
            0)] <- x

        colnames(cor) <- situation.states.time
        rownames(cor) <- situation.states.time

        return(cor)
    }


    ## La fonction recode() fournit la décomposition en composantes principales de
    ## l'espace des profils Les coordonnées de chaque situation sur les axes
    ## principaux La décomposition de chaque sequence sur les axes principaux La
    ## distance entre les séquences peut alors être définie comme la distance
    ## euclidienne

    recode <- function(disj, prin) {
        recode <- disj %*% prin
        return(recode)
    }


    ## Calculs

    emlt <- list()
    emlt$states <- alphabet(seqdata)
    emlt$period <- dim(seqdata)[2]
    emlt$sit.time <- rep(1:dim(seqdata)[2], each = length(alphabet(seqdata)))
    emlt$situations <- paste(alphabet(seqdata), emlt$sit.time, sep = "")
    emlt$sit.states <- rep(alphabet(seqdata), times = dim(seqdata)[2])
    emlt$sit.freq <- freq(seqdata)

    emlt$a <- disjonctif(seqdata)
    if (!weighted) {
        attr(seqdata, "weights") <- NULL
    }
    emlt$sit.transrate <- transrate(seqdata, emlt$a, poids = attr(seqdata, "weights"))
    emlt$sit.profil <- profil(seqdata, emlt$sit.transrate, param_a = a, param_b = b)
    emlt$c <- distsquare(seqdata, emlt$sit.profil)
    emlt$d <- benz(seqdata, emlt$c)
    emlt$pca <- princomp(emlt$d, cor = TRUE)
    emlt$coord <- recode(emlt$a[, which(emlt$sit.freq != 0)], emlt$pca$scores)
    emlt$sit.cor <- corel(seqdata, emlt$d)
    emlt$seqdata <- seqdata
    class(emlt) <- "emlt"

    return(emlt)
}

print.emlt <- function(x, ...) {
    ## Print some useful things
    cat("Periods of observation: ", x$period, "\n")
    cat("State alphabet: ")
    cat(x$state,"\n\n")
    cat("Frequencies of situations (time-stamped states):  \n")
    print(x$sit.freq)
    cat("\n")
    cat("sit.cor: Correlations between situations\n         (proximities/links between situations) \n\n")
    cat("coord: Euclidean coordinates of the sequences\n       (can be used as input to clustering algorithms)")
}

plot.emlt <- function(x, from, to, delay = NULL, leg = TRUE, type = "cor", cex = 0.7,
    compx = 1, compy = 2, ...) {
    if (type == "cor") {
        delai <- 0
        subtitle <- NULL
        if (!is.null(delay)) {
            delai <- delay
            subtitle <- paste("with a delay=", delay)
        }

        period <- x$period
        par(new = F)
        axe1d <- 0
        axe1by <- floor((period - delai + 16)/4)
        axe1f <- axe1d + 4 * (axe1by)
        listto <- which(x$sit.states == to)[delai:period]
        listfrom <- which(x$sit.states == from[1])[1:(period - delai)]
        title1 <- paste(from, collapse = ",")
        title <- c(paste("correlation between being in ", title1, "and being in ",
            to), subtitle)
        min <- min(diag(x$sit.cor[listfrom, listto]), na.rm = TRUE)
        max <- max(diag(x$sit.cor[listfrom, listto]), na.rm = TRUE)
        for (i in 1:length(from)) {
            listfrom <- which(x$sit.states == from[i])[1:(period - delai)]

            mi <- min(diag(x$sit.cor[listfrom, listto]), na.rm = TRUE)

            ma <- max(diag(x$sit.cor[listfrom, listto]), na.rm = TRUE)

            if (mi < min) {
                min <- mi
            }
            if (ma > max) {
                max <- ma
            }

        }


        for (i in 1:length(from)) {
            listfrom <- which(x$sit.states == from[i])[1:(period - delai)]
            par(lty = i)
            if (i == 1) {
                plot(diag(x$sit.cor[listfrom, listto]), type = "l", main = title,
                  xlim = c(0, period), ylim = c(min, max), xlab = "time", ylab = "correlation",
                  xaxp = c(axe1d, axe1f, by = axe1by))
            }
            if (i > 1) {
                points(diag(x$sit.cor[listfrom, listto]), type = "l", xlab = "time",
                  ylab = "correlation", xaxp = c(axe1d, axe1f, by = axe1by), xlim = c(0,
                    period), ylim = c(mi, max))
            }


        }
        if (leg) {
            legend(x = c(period/3, period/2), y = c(min + 0.7 * (max - min), max),
                bty = "n", legend = from, lty = 1:length(from))
        }


    }
    if (type == "pca") {
        title <- c(paste("PCA - Principal plane of situations space - components",
            compx, "and ", compy))
        ptext <- x$situations[which(x$sit.freq != 0)]
        plot(x$pca$scores[, compx], x$pca$scores[, compy], xlab = paste("comp", compx),
            ylab = paste("comp", compy), main = title)
        text(x$pca$scores[, compx], x$pca$scores[, compy], ptext, cex = cex)
    }

}
