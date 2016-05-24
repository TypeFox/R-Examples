concordance<-function (db, test = "Default", B = 1000, alpha = 0.05) 
{
    db2 <- db * (db - 1)
    sumcol <- colSums(db)
    sumrow <- rowSums(db)
    tot <- sum(db)
    vet <- list()
    pij <- matrix(, nrow = nrow(db), ncol = ncol(db))
    for (i in 1:length(sumrow)) {
        for (j in 1:length(sumcol)) {
            pij[i, j] <- db2[i, j]/(sumrow[i] * (sumrow[i] - 
                1))
        }
    }
    pi <- rowSums(pij)
    p <- mean(pi)
    pj <- sumcol/tot
    pj2 <- pj^2
    pe <- sum(pj2)
    fleiss.kappa <- (p - pe)/(1 - pe)
    s <- (ncol(db) * p - 1)/(ncol(db) - 1)
    s.boot <- c()
    pi.v.boot <- replicate(B, pi.boot <- sample(pi, size = nrow(db), 
        replace = TRUE))
    p.boot <- apply(pi.v.boot, 2, mean)
    for (i in 1:B) {
        s.boot[i] <- (ncol(db) * p.boot[i] - 1)/(ncol(db) - 1)
    }
    s.boot.ci <- quantile(s.boot, probs = c(alpha/2, 1 - alpha/2))
    Default <- function(vet) {
        s.vet <- c(s, s.boot.ci)
        names(s.vet) <- c("S", "LCL", "UCL")
        vet <- list(Fleiss = fleiss.kappa, Statistic = s.vet)
        cat(paste("Inter-rater Agreement"), "\n")
        vet
    }
    Normal <- function(db) {
        stat.test <- s * nrow(db) * sqrt((ncol(db) - 1)/(2 * 
            sum(1/(rowSums(db) * (rowSums(db) - 1)))))
        bin <- stat.test > qnorm(1 - alpha)
        pvalue <- 1 - pnorm(stat.test)
        s.vet <- c(s, s.boot.ci, pvalue)
        names(s.vet) <- c("S", "LCL", "UCL", "pvalue")
        vet <- list(Fleiss = fleiss.kappa, Statistic = s.vet)
        cat(paste("Inter-rater Agreement"), "\n")
        vet
    }
    MC <- function(db) {
        matrix.mc <- replicate(B, matrix.rep <- matrix(, nrow = nrow(db), 
            ncol = ncol(db)))
        pij.mc <- replicate(B, pij.rep <- matrix(, nrow = nrow(db), 
            ncol = ncol(db)))
        for (k in 1:B) {
            for (j in 1:nrow(db)) matrix.mc[j, , k] <- t(rmultinom(1, 
                size = rowSums(db)[j], prob = rep(1/(ncol(db)), 
                  ncol(db))))
        }
        for (k in 1:B) {
            for (i in 1:nrow(db)) {
                for (j in 1:ncol(db)) {
                  pij.mc[i, j, k] <- matrix.mc[i, j, k] * (matrix.mc[i, 
                    j, k] - 1)/(rowSums(db)[i] * (rowSums(db)[i] - 
                    1))
                }
            }
        }
        pi.mc <- list()
        for (i in 1:B) {
            pi.mc[[i]] <- apply(pij.mc[, , i], 1, sum)
        }
        p.mc <- c()
        for (k in 1:B) {
            p.mc[k] <- mean(pi.mc[[k]])
        }
        s.mc <- c()
        for (k in 1:B) {
            s.mc[k] <- ((p.mc[k] * ncol(db)) - 1)/(ncol(db) - 
                1)
        }
        crit.s.mc <- quantile(s.mc, 0.95)
        binary <- c()
        for (i in 1:length(s.mc)) {
            binary[i] <- (s.mc[i] >= s)
        }
        pvalue <- sum(binary)/B
        s.vet <- c(s, s.boot.ci, pvalue)
        names(s.vet) <- c("S", "LCL", "UCL", "pvalue")
        vet <- list(Fleiss = fleiss.kappa, Statistic = s.vet)
        cat(paste("Inter-rater Agreement"), "\n")
        vet
    }
    if (length(unique(rowSums(db))) == 1) {
        std.kappaj <- sqrt(2/(nrow(db) * (rowSums(db)[1] * (rowSums(db)[1] - 
            1))))
		pj.1<-1-pj
		pj.2<-1-2*pj
		std.num<-sqrt((sum(pj*pj.1))^2-sum(pj*pj.1*pj.2))
		std.den<-sum(pj*pj.1)
		std.adj<-std.num/std.den
		std.kappa<-std.kappaj*std.adj
        ci.low <- fleiss.kappa - qnorm(1 - alpha/2) * std.kappa
        ci.upp <- fleiss.kappa + qnorm(1 - alpha/2) * std.kappa
        z.fleiss <- fleiss.kappa/std.kappa
        pvalue.fleiss <- 1 - pnorm(z.fleiss)
        if (ci.upp > 1 | ci.low < (-1/(rowSums(db)[1] - 1))) {
            ci.low <- (1/(rowSums(db)[1] - 1))
            ci.upp <- 1
        }
        Chisq <- function(db) {
            stat.test <- nrow(db) * (ncol(db) - 1) * ((rowSums(db) - 
                1)[1] * s + 1)
            bin <- stat.test > qchisq(alpha, df = (nrow(db) * 
                (ncol(db) - 1)))
            pvalue <- 1 - pchisq(stat.test, df = (nrow(db) * 
                (ncol(db) - 1)))
            kappa.vet <- c(fleiss.kappa, ci.low, ci.upp, std.kappa, 
                z.fleiss, pvalue.fleiss)
            names(kappa.vet) <- c("Kappa", "LCL", "UCL", "Std.Error", 
                "Z value", "Pr(>|z|)")
            s.vet <- c(s, s.boot.ci, pvalue)
            names(s.vet) <- c("S", "LCL", "UCL", "pvalue")
            vet <- list(Fleiss = kappa.vet, Statistic = s.vet)
            cat(paste("Inter-rater Agreement"), "\n")
            vet
        }
        Normal2 <- function(db) {
            stat.test <- s * nrow(db) * sqrt((ncol(db) - 1)/(2 * 
                sum(1/(rowSums(db) * (rowSums(db) - 1)))))
            bin <- stat.test > qnorm(1 - alpha)
            pvalue <- 1 - pnorm(stat.test)
            kappa.vet <- c(fleiss.kappa, ci.low, ci.upp, std.kappa, 
                z.fleiss, pvalue.fleiss)
            names(kappa.vet) <- c("Kappa", "LCL", "UCL", "Std.Error", 
                "Z value", "Pr(>|z|)")
            s.vet <- c(s, s.boot.ci, pvalue)
            names(s.vet) <- c("S", "LCL", "UCL", "pvalue")
            vet <- list(Fleiss = kappa.vet, Statistic = s.vet)
            cat(paste("Inter-rater Agreement"), "\n")
            vet
        }
        MC2 <- function(db) {
            matrix.mc <- replicate(B, matrix.rep <- t(rmultinom(nrow(db), 
                size = rowSums(db)[1], prob = rep(1/(ncol(db)), 
                  ncol(db)))))
            pij.mc <- replicate(B, pij.rep <- matrix(, nrow = nrow(db), 
                ncol = ncol(db)))
            for (k in 1:B) {
                for (i in 1:nrow(db)) {
                  for (j in 1:ncol(db)) {
                    pij.mc[i, j, k] <- matrix.mc[i, j, k] * (matrix.mc[i, 
                      j, k] - 1)/(rowSums(db)[i] * (rowSums(db)[i] - 
                      1))
                  }
                }
            }
            pi.mc <- list()
            for (i in 1:B) {
                pi.mc[[i]] <- apply(pij.mc[, , i], 1, sum)
            }
            p.mc <- c()
            for (k in 1:B) {
                p.mc[k] <- mean(pi.mc[[k]])
            }
            s.mc <- c()
            for (k in 1:B) {
                s.mc[k] <- ((p.mc[k] * ncol(db)) - 1)/(ncol(db) - 
                  1)
            }
            crit.s.mc <- quantile(s.mc, 0.95)
            binary <- c()
            for (i in 1:length(s.mc)) {
                binary[i] <- (s.mc[i] >= s)
            }
            pvalue <- sum(binary)/B
            kappa.vet <- c(fleiss.kappa, ci.low, ci.upp, std.kappa, 
                z.fleiss, pvalue.fleiss)
            names(kappa.vet) <- c("Kappa", "LCL", "UCL", "Std.Error", 
                "Z value", "Pr(>|z|)")
            s.vet <- c(s, s.boot.ci, pvalue)
            names(s.vet) <- c("S", "LCL", "UCL", "pvalue")
            vet <- list(Fleiss = kappa.vet, Statistic = s.vet)
            cat(paste("Inter-rater Agreement"), "\n")
            vet
        }
        switch(test, Normal = Normal2(db), MC = MC2(db), Chisq = Chisq(db), 
            Default = Default(vet))
    }
    else {
        switch(test, Normal = Normal(db), MC = MC(db), Default = Default(vet))
    }
}