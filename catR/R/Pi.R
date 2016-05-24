Pi<-function (th, it, model = NULL, D = 1) 
{
    it <- rbind(it)
    if (is.null(model)) {
        a <- it[, 1]
        b <- it[, 2]
        c <- it[, 3]
        d <- it[, 4]
        e <- exp(D * a * (th - b))
        Pi <- c + (d - c) * e/(1 + e)
        Pi[Pi == 0] <- 1e-10
        Pi[Pi == 1] <- 1 - 1e-10
        dPi <- D * a * e * (d - c)/(1 + e)^2
        d2Pi <- D^2 * a^2 * e * (1 - e) * (d - c)/(1 + e)^3
        d3Pi <- D^3 * a^3 * e * (d - c) * (e^2 - 4 * e + 1)/(1 + 
            e)^4
        res <- list(Pi = Pi, dPi = dPi, d2Pi = d2Pi, d3Pi = d3Pi)
    }
    else {
        if (sum(model == c("GRM", "MGRM", "PCM", "GPCM", "RSM", 
            "NRM")) == 0) 
            stop("invalid 'model' name'", call. = FALSE)
        if (model == "GRM" | model == "MGRM") {
            if (model == "GRM") 
                prov <- prov1 <- prov2 <- prov3 <- matrix(NA, 
                  nrow(it), ncol(it))
            else prov <- prov1 <- prov2 <- prov3 <- matrix(NA, 
                nrow(it), ncol(it) - 1)
            for (i in 1:nrow(it)) {
                aj <- it[i, 1]
                if (model == "GRM") 
                  bj <- it[i, 2:ncol(it)]
                else bj <- it[i, 2] - it[i, 3:ncol(it)]
                bj <- bj[!is.na(bj)]
                ej <- exp(D * aj * (th - bj))
                Pjs <- ej/(1 + ej)
                Pjs <- c(1, Pjs, 0)
                dPjs <- D * aj * Pjs * (1 - Pjs)
                d2Pjs <- D * aj * (dPjs - 2 * Pjs * dPjs)
                d3Pjs <- D * aj * (d2Pjs - 2 * dPjs^2 - 2 * Pjs * 
                  d2Pjs)
                n <- length(Pjs)
                prov[i, 1:(n - 1)] <- Pjs[1:(n - 1)] - Pjs[2:n]
                prov1[i, 1:(n - 1)] <- dPjs[1:(n - 1)] - dPjs[2:n]
                prov2[i, 1:(n - 1)] <- d2Pjs[1:(n - 1)] - d2Pjs[2:n]
                prov3[i, 1:(n - 1)] <- d3Pjs[1:(n - 1)] - d3Pjs[2:n]
            }
        }
        else {
            nc <- switch(model, PCM = ncol(it) + 1, GPCM = ncol(it), 
                RSM = ncol(it), NRM = ncol(it)/2 + 1)
            prov <- prov1 <- prov2 <- prov3 <- matrix(NA, nrow(it), 
                nc)
            for (i in 1:nrow(it)) {
                dj <- v<-0
                if (model == "PCM") {
                  for (t in 1:ncol(it)){
 dj <- c(dj, dj[t] + D * (th - it[i, t]))
v<-c(v,t)
}
                }
                if (model == "GPCM") {
                  for (t in 1:(ncol(it) - 1)){
 dj <- c(dj, dj[t] + it[i, 1] * D * (th - it[i, t + 1]))
v<-c(v,it[i, 1]*t)
}
                }
                if (model == "RSM") {
                  for (t in 1:(ncol(it) - 1)){
 dj <- c(dj, dj[t] + D * (th - (it[i, 1] + it[i, t + 1])))
v<-c(v,t)
}
                }
                if (model == "NRM") {
                  for (t in 1:(ncol(it)/2)){
 dj <- c(dj, it[i, (2 * (t - 1) + 1)] * th + it[i, (2 * t)])
v<-c(v,it[i, (2 * (t - 1) + 1)])
}
                }
v<-v[!is.na(dj)]
                dj <- dj[!is.na(dj)]
                Gammaj <- exp(dj)
                dGammaj <- Gammaj * v
                d2Gammaj <- Gammaj * v^2
                d3Gammaj <- Gammaj * v^3
                Sg <- sum(Gammaj)
                Sdg <- sum(dGammaj)
                Sd2g <- sum(d2Gammaj)
                Sd3g <- sum(d3Gammaj)
                n <- length(Gammaj)
                prov[i, 1:n] <- Gammaj/Sg
                prov1[i, 1:n] <- dGammaj/Sg - Gammaj * Sdg/Sg^2
                prov2[i, 1:n] <- d2Gammaj/Sg - 2 * dGammaj * 
                  Sdg/Sg^2 - Gammaj * Sd2g/Sg^2 + 2 * Gammaj * 
                  Sdg^2/Sg^3
                prov3[i, 1:n] <- d3Gammaj/Sg - (Gammaj * Sd3g + 
                  3 * dGammaj * Sd2g + 3 * d2Gammaj * Sdg)/Sg^2 + 
                  (6 * Gammaj * Sdg * Sd2g + 6 * dGammaj * Sdg^2)/Sg^3 - 
                  6 * Gammaj * Sdg^3/Sg^4
            }
        }
        cn <- "cat0"
        for (i in 1:(ncol(prov) - 1)) cn <- c(cn, paste("cat", 
            i, sep = ""))
        colnames(prov) <- colnames(prov1) <- colnames(prov2) <- colnames(prov3) <- cn
        rn <- "Item1"
        if (nrow(it) > 1) {
            for (i in 2:nrow(prov)) rn <- c(rn, paste("Item", 
                i, sep = ""))
        }
        rownames(prov) <- rownames(prov1) <- rownames(prov2) <- rownames(prov3) <- rn
        res <- list(Pi = prov, dPi = prov1, d2Pi = prov2, d3Pi = prov3)
    }
    return(res)
}


