`AMMI` <-
function(ENV,GEN,REP,Y,MSE = 0,console=FALSE,PC=FALSE)
{
    name.y <- paste(deparse(substitute(Y)))
    if(console)cat("\nANALYSIS AMMI: ", name.y, "\nClass level information\n")
    ENV <- as.factor(ENV)
    GEN <- as.factor(GEN)
    nenv <- length(unique(ENV))
    ngen <- length(unique(GEN))
    if(console)cat("\nENV: ", unique(as.character(ENV)))
    if(console)cat("\nGEN: ", unique(as.character(GEN)))
    minimo <- min(ngen, nenv)
    if (length(REP) > 1) {
        REP <- as.factor(REP)
        nrep <- length(unique(REP))
        if(console)cat("\nREP: ", unique(REP))
        if(console)cat("\n\nNumber of observations: ", length(na.omit(Y)),
            "\n\n")
        modelo <- aov(Y ~ ENV + REP %in% ENV + GEN + ENV:GEN)
        if(console)cat("model Y:", name.y, " ~ ENV + REP%in%ENV + GEN + ENV:GEN\n")
        if(console)cat("Random effect REP%in%ENV\n\n")
        mm <- anova(modelo)
        nn <- mm[2, ]
        mm[2, ] <- mm[3, ]
        mm[3, ] <- nn
        row.names(mm)[2] <- "REP(ENV)"
        row.names(mm)[3] <- "GEN     "
        mm[1, 4] <- mm[1, 3]/mm[2, 3]
        mm[1, 5] <- 1 - pf(mm[1, 4], mm[1, 1], mm[2, 1])
        if(console)print(mm)
        anova <- mm
        DFE <- df.residual(modelo)
        MSE <- deviance(modelo)/DFE
        medy <- mean(Y, na.rm = TRUE)
        if(console)cat("\nCoeff var", "\tMean", name.y, "\n")
        if(console)cat(sqrt(MSE) * 100/medy, "\t", medy, "\n")
    }
    else {
        DFE <- nenv * (ngen - 1) * (REP - 1)
        DFEa <- nenv * (REP - 1)
        nrep <- REP
        modelo <- aov(Y ~ ENV + GEN)
        xx <- as.matrix(anova(modelo))
        xx <- rbind(xx[1, ], xx[1, ], xx[2:3, ],xx[3,])
		row.names(xx)[4] <- "ENV:GEN"
		row.names(xx)[5] <- "Residuals"
        xx[2, 1] <- DFEa
        xx[2, 2:5] <- NA
        xx[, 2] <- xx[, 2] * nrep
        xx[, 3] <- xx[, 3] * nrep
        xx[5, 1] <- DFE
        xx[5, 3] <- MSE
        xx[5, 2] <- MSE * DFE
        xx[1, 4] <- NA
        if(MSE>0) xx[3, 4] <- xx[3, 3]/MSE
        if(MSE>0)xx[4, 4] <- xx[4, 3]/MSE
        xx[1, 5] <- NA
        if(DFE >0){
        xx[3, 5] <- 1 - pf(xx[3, 4], xx[3, 1], DFE)
        xx[4, 5] <- 1 - pf(xx[4, 4], xx[4, 1], DFE)
        }
        row.names(xx)[1] <- "ENV     "
        row.names(xx)[2] <- "REP(ENV)"
        if(console)cat("\nREP: ", REP)
        if(console)cat("\n\nNumber of means: ", length(na.omit(Y)), "\n")
        if(console)cat("\nDependent Variable:", name.y, "\n\nAnalysis of variance\n")
        mm<-xx
        if(console)print(xx, na.print = "")
        medy <- mean(Y, na.rm = TRUE)
        if(console)cat("\nCoeff var", "\tMean", name.y, "\n")
        if(console)cat(sqrt(MSE) * 100/medy, "\t", medy, "\n")
    }

    raw <- data.frame(ENV, GEN, Y)
    MEDIAS <- tapply(raw[, 3], raw[, c(1, 2)], mean)
    xx <- rownames(MEDIAS)
    yy <- colnames(MEDIAS)
    fila <- length(xx)
    col <- length(yy)
    total <- fila * col
    x <- character(length = total)
    y <- character(length = total)
    z <- numeric(length = total)
    k <- 0
    for (i in 1:fila) {
        for (j in 1:col) {
            k <- k + 1
            x[k] <- xx[i]
            y[k] <- yy[j]
            z[k] <- MEDIAS[i, j]
        }
    }
    MEDIAS <- data.frame(ENV = x, GEN = y, Y = z)
    x <- MEDIAS[, 1]
    y <- MEDIAS[, 2]
    z <- MEDIAS[, 3]
    modelo2 <- lm(z ~ x + y)
    for (i in 1:length(z)) {
        if (is.na(z[i]))
            z[i] <- predict(modelo2, data.frame(x = MEDIAS[i,
                1], y = MEDIAS[i, 2]))
    }
    MEDIAS <- data.frame(ENV = x, GEN = y, Y = z)
    modelo1 <- lm(Y ~ ENV + GEN, data = MEDIAS)
    residual <- modelo1$residuals
    MEDIAS <- data.frame(MEDIAS, RESIDUAL = residual)
    mlabel <- names(MEDIAS)
    names(MEDIAS) <- c(mlabel[1:2], name.y, mlabel[4])
    OUTRES <- MEDIAS[order(MEDIAS[, 1], MEDIAS[, 2]), ]
    OUTRES2 <- by(OUTRES[, 4], OUTRES[, c(2, 1)], function(x) sum(x,na.rm=TRUE))
    OUTMED <- by(OUTRES[, 3], OUTRES[, c(2, 1)], function(x) sum(x,na.rm=TRUE))
    s <- svd(OUTRES2)
    U <- s$u
    L <- s$d
    V <- s$v
    L <- L[1:minimo]
    SS <- (L^2) * nrep
    SUMA <- sum(SS)
    percent <- (1/SUMA)*SS*100
#    minimo <- min(ngen, nenv)
    DFAMMI <- rep(0, minimo)
    acum <- DFAMMI
    MSAMMI <- DFAMMI
    F.AMMI <- DFAMMI
    PROBF <- DFAMMI
    acumula <- 0
    for (i in 1:(minimo-1)) {
        DF <- (ngen - 1) + (nenv - 1) - (2 * i - 1)
        if (DF <= 0)
            break
        DFAMMI[i] <- DF
        acumula <- acumula + percent[i]
        acum[i] <- acum[i] + acumula
        MSAMMI[i] <- SS[i]/DFAMMI[i]
        if(MSE>0)F.AMMI[i] <- round(MSAMMI[i]/MSE, 2)
        else F.AMMI[i] <-NA
        if(DFE>0)
        PROBF[i] <- round(1 - pf(F.AMMI[i], DFAMMI[i], DFE), 4)
        else PROBF[i]<-NA
    }
	percent<-round(percent,1)
	acum<-round(acum,1)
    SS <- round(SS, 6)
    MSAMMI <- round(MSAMMI, 6)
    SSAMMI <- data.frame(percent, acum, Df = DFAMMI, `Sum Sq` = SS,
        `Mean Sq` = MSAMMI, `F value` = F.AMMI, Pr.F = PROBF)
    nssammi <- nrow(SSAMMI)
    SSAMMI <- SSAMMI[SSAMMI$Df > 0, ]
    nss <- nrow(SSAMMI)
    row.names(SSAMMI) <- paste("PC", 1:nss, sep = "")
    if(console){cat("\nAnalysis\n")
    print(SSAMMI)
    }
    LL <- sqrt(diag(L))
    SCOREG <- U %*% LL
    SCOREE <- V %*% LL
    SCORES <- rbind(SCOREG, SCOREE)
    colnames(SCORES) <- paste("PC", 1:nssammi, sep = "")
    MSCORES <- SCORES[1:ngen, ]
    NSCORES <- SCORES[(ngen + 1):(ngen + nenv), ]
    MGEN <- data.frame(type = "GEN", Y = apply(OUTMED, 1, mean),
        MSCORES)
    MENV <- data.frame(type = "ENV", Y = apply(OUTMED, 2, mean),
        NSCORES)
    bplot <- rbind(MGEN, MENV)
    bplot <- bplot[, 1:(nss + 2)]
    mlabel <- names(bplot)
    names(bplot) <- c(mlabel[1], name.y, mlabel[c(-1, -2)])
if( minimo <= 2) {
    cat("\nWarning. The analysis AMMI is not possible.")
    cat("\nThe number of environments and number of genotypes must be greater than 2\n")
    }
if(PC) PC<- princomp(OUTRES2, cor = FALSE)
object<-list(ANOVA=mm,genXenv = OUTRES2, analysis = SSAMMI, means = MEDIAS, biplot = bplot,PC=PC)
class(object)<-"AMMI"
invisible(object)
}
