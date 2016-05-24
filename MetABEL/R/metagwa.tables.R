##' Pairwise meta-analysis of results from GWA scans
##'
##' Performes meta-analysis of results of two individual GWA studies
##'
##' Original data frames should contain a number of variables,
##' such as allelic coding, code for the effective allele,
##' etc. Please refer to the \code{formetascore} function of the
##' \code{GenABEL} package for details.
##' @param data.x First data frame with GWA data
##' @param data.y Second data frame with GWA data
##' @param name.x First study name
##' @param name.y Second study name
##' @param precorrect Should GC be applied to the original data before pooling
##' @param correct.pooled Whether to apply Genomic Control correction to
##' the study named "POOLED"
##' @return   A data frame containing summary statistics and essential details
##'  of the individual studies
##' @author Yurii Aulchenko
##' @seealso \code{\link{metagwa.files}}
##' @keywords htest
##' @export
"metagwa.tables" <-
    function(data.x, data.y, name.x="P1", name.y="P2",
             precorrect=TRUE, correct.pooled=FALSE) {

        required <- c("name", "allele1", "allele2", "effallele", "beta", "sebeta")
        optional <- c("chromosome", "position", "n", "effallelefreq")
        n.x <- names(data.x)
        essen <- match(required, n.x)

        if (any(is.na(essen))) {
            cat("Required fields missing in data.x:",
                required[is.na(essen)], "\n")
            stop()
        }

        opt <- match(optional, n.x)

        if (any(is.na(essen))) {
            cat("Optional fields missing in data.x:",
                required[is.na(essen)], "\n")
        }

        n.y <- names(data.y)
        essen <- match(required, n.y)

        if (any(is.na(essen))) {
            cat("Required fields missing in data.y:",
                required[is.na(essen)], "\n");
            stop()
        }

        opt <- match(optional, n.y)
        if (any(is.na(essen))) {
            cat("Optional fields missing in data.y:",
                required[is.na(essen)], "\n");
        }

        pos.x <- 1
        if (is.na(match("position", names(data.x)))) pos.x <- 0
        pos.y <- 1
        if (is.na(match("position", names(data.y)))) pos.y <- 0
        chr.x <- 1
        if (is.na(match("chromosome", names(data.x)))) chr.x <- 0
        chr.y <- 1
        if (is.na(match("chromosome", names(data.y)))) chr.y <- 0

        df <- merge(data.x, data.y, by="name", all.x=TRUE, all.y=TRUE)

        ## Do checks
        err1 <- which(df$allele1.x != df$effallele.x &
                      df$allele2.x != df$effallele.x)
        if (length(err1) > 0) {
                cat("Effective allele does not correspond to either allele1 or allele2 in",
                    name.x, "for SNPs\n")
                print(df$name[err1])
                stop("Stopped")
        }

        err2 <- which(df$allele1.y != df$effallele.y &
                      df$allele2.y != df$effallele.y)
        if (length(err1) > 0) {
                cat("Effective allele does not correspond to either allele1 or allele2 in",
                    name.x, "for SNPs\n")
                print(df$name[err1])
                stop("Stopped")
        }
        gc()

        ## OK, now check that effective allele is the second one
        swall1 <- which(df$allele2.x != df$effallele.x)
        if (length(swall1) > 0) {
                keep <- df$allele2.x[swall1]
                df$allele2.x[swall1] <- df$allele1.x[swall1]
                df$allele1.x[swall1] <- keep
        }

        swall2 <- which(df$allele2.y != df$effallele.y)
        if (length(swall2) > 0) {
                keep <- df$allele2.y[swall2]
                df$allele2.y[swall2] <- df$allele1.y[swall2]
                df$allele1.y[swall2] <- keep
        }
        gc()

        ## OK, now put TA->AT, GA->AG etc (A>T>G>C)
        recX <- which(df$allele1.x=="T" & df$allele2.x=="A")
        if (length(recX)>0) df <- swap.x(df, recX)
        recX <- which(df$allele1.x=="G" & df$allele2.x=="A")
        if (length(recX)>0) df <- swap.x(df, recX)
        recX <- which(df$allele1.x=="C" & df$allele2.x=="A")
        if (length(recX)>0) df <- swap.x(df, recX)
        recX <- which(df$allele1.x=="G" & df$allele2.x=="T")
        if (length(recX)>0) df <- swap.x(df, recX)
        recX <- which(df$allele1.x=="C" & df$allele2.x=="T")
        if (length(recX)>0) df <- swap.x(df, recX)
        recX <- which(df$allele1.x=="C" & df$allele2.x=="G")
        if (length(recX)>0) df <- swap.x(df, recX)
        recX <- which(df$allele1.y=="T" & df$allele2.y=="A")
        if (length(recX)>0) df <- swap.y(df, recX)
        recX <- which(df$allele1.y=="G" & df$allele2.y=="A")
        if (length(recX)>0) df <- swap.y(df, recX)
        recX <- which(df$allele1.y=="C" & df$allele2.y=="A")
        if (length(recX)>0) df <- swap.y(df, recX)
        recX <- which(df$allele1.y=="G" & df$allele2.y=="T")
        if (length(recX)>0) df <- swap.y(df, recX)
        recX <- which(df$allele1.y=="C" & df$allele2.y=="T")
        if (length(recX)>0) df <- swap.y(df, recX)
        recX <- which(df$allele1.y=="C" & df$allele2.y=="G")
        if (length(recX)>0) df <- swap.y(df, recX)
        gc()

        ## OK, now try to unify coding between studies
        swcode <- which(df$effallele.x != df$effallele.y |
                        df$allele1.x != df$allele1.y |
                        df$allele2.x != df$allele2.y)
        if (length(swcode) > 0) {
                df <- recode.y(df, swcode)
                recX <- which(df$allele1.y=="T" & df$allele2.y=="A")
                if (length(recX)>0) df <- swap.y(df, recX)
                recX <- which(df$allele1.y=="G" & df$allele2.y=="A")
                if (length(recX)>0) df <- swap.y(df, recX)
                recX <- which(df$allele1.y=="C" & df$allele2.y=="A")
                if (length(recX)>0) df <- swap.y(df, recX)
                recX <- which(df$allele1.y=="G" & df$allele2.y=="T")
                if (length(recX)>0) df <- swap.y(df, recX)
                recX <- which(df$allele1.y=="C" & df$allele2.y=="T")
                if (length(recX)>0) df <- swap.y(df, recX)
                recX <- which(df$allele1.y=="C" & df$allele2.y=="G")
                if (length(recX)>0) df <- swap.y(df, recX)
        }
        gc()

        ## OK, now nullify ambigous coding with wrong strand info
        toerX <- which(df$allele1.x == "A" & df$allele2.x == "T" &
                       (df$strand.x!="+" & df$strand.x!="-" &
                        df$strand.x!="TOP" & df$strand.x!="BOT"))
        if (length(toerX) > 0) {
                cat("AT coding without +/-/TOP/BOT strand info in",
                    name.x, "\n")
                cat(length(toerX), "SNPs nullified\n")
                df$n.x[toerX] <- NA
                df$beta.x[toerX] <- NA
        }
        gc()

        toerX <- which(df$allele1.x == "G" & df$allele2.x == "C" &
                       (df$strand.x!="+" & df$strand.x!="-" &
                        df$strand.x!="TOP" & df$strand.x!="BOT"))
        if (length(toerX) > 0) {
                cat("GC coding without +/-/TOP/BOT strand info in",
                    name.x, "\n")
                cat(length(toerX), "SNPs nullified\n")
                df$n.x[toerX] <- NA
                df$beta.x[toerX] <- NA
        }
        gc()

        toerX <- which(df$allele1.y == "A" & df$allele2.y == "T" &
                       (df$strand.y!="+" & df$strand.y!="-" &
                        df$strand.y!="TOP" & df$strand.y!="BOT"))
        if (length(toerX) > 0) {
                cat("AT coding without +/-/TOP/BOT strand info in",
                    name.y, "\n")
                cat(length(toerX), "SNPs nullified\n")
                df$n.y[toerX] <- NA
                df$beta.y[toerX] <- NA
        }
        gc()

        toerX <- which(df$allele1.y == "G" & df$allele2.y == "C" &
                       (df$strand.y!="+" & df$strand.y!="-" &
                        df$strand.y!="TOP" & df$strand.y!="BOT"))
        if (length(toerX) > 0) {
                cat("GC coding without +/-/TOP/BOT strand info in",
                    name.y, "\n")
                cat(length(toerX), "SNPs nullified\n")
                df$n.y[toerX] <- NA
                df$beta.y[toerX] <- NA
        }
        gc()

        swcode <- which(df$effallele.x != df$effallele.y |
                        df$allele1.x != df$allele1.y |
                        df$allele2.x != df$allele2.y)
        if (length(swcode) > 0) {
                cat("Different coding in two populations\n")
                cat(length(swcode), "SNPs removed\n")
                tmp <- df[swcode, ]
                write.csv(tmp, file="failedcode.csv", row.names=FALSE)
                df <- df[!(c(1:(dim(df)[1])) %in% swcode), ]
        }
        gc()

        swcode <- which(is.na(df$beta.x) & is.na(df$beta.y))
        if (length(swcode)>0) {
                cat("NA for betas in both populations\n")
                cat(length(swcode), "SNPs removed\n")
                df <- df[!(c(1:(dim(df)[1])) %in% swcode), ]
        }
        gc()

        ## OK
        cat("analysing ... \n")
#	if (pop==pops[2]) name.x <- pops[1] else name.x <- "POOLED"
        metaout <- dometa(df, name.x=name.x, name.y=name.y,
                          precorrect=precorrect,
                          correct.pooled=correct.pooled,
                          pos.x, pos.y, chr.x, chr.y)
        cat("... DONE\n")
        gc()
        metaout
}


dometa <- function(data, name.x, name.y, precorrect=FALSE,
                   correct.pooled=FALSE, pos.x, pos.y, chr.x, chr.y) {
#	prop <- 0.98
        chi2med <- qchisq(.5, 1)
        b.x <- (data$beta.x)
        b.y <- (data$beta.y)
        se.x <- abs(data$sebeta.x)
        se.y <- abs(data$sebeta.y)

        if (name.x == "POOLED" & !correct.pooled) {
            lam.x <- 1.0
        } else {
            lam.x <- median(b.x * b.x / (se.x * se.x), na.rm=TRUE) / chi2med
        }

        lam.y <- median(b.y * b.y / (se.y * se.y), na.rm=TRUE) / chi2med

        if (name.x != "POOLED")
            cat("Lambda", name.x, "=", lam.x, "\n")
        else if (correct.pooled)
            cat("Lambda POOLED previous =", lam.x, "\n")

        cat("Lambda", name.y, "=", lam.y, "\n")
        if (lam.x < 1) {
            warning(paste("Lambda", name.x, "< 1; constrained to 1"),
                    immediate.=TRUE); lam.x <- 1
        }
        if (lam.y < 1) {
            warning(paste("Lambda", name.y, "< 1; constrained to 1"),
                    immediate.=TRUE); lam.y <- 1
        }
        if (precorrect & name.x!="POOLED") {
            se.x <- sqrt(se.x*se.x*lam.x)
        }
        if (precorrect) {
            se.y <- sqrt(se.y*se.y*lam.y)
        }

        lam.x <- median(b.x * b.x / (se.x * se.x), na.rm=TRUE)/chi2med
        lam.y <- median(b.y * b.y / (se.y * se.y), na.rm=TRUE)/chi2med

        if (name.x != "POOLED")
            cat("Corrected Lambda", name.x, "=", lam.x, "\n")
        else if (correct.pooled)
            cat("Corrected Lambda POOLED previous =", lam.x, "\n")

        cat("Corrected Lambda", name.y, "=", lam.y, "\n")
        w2.x     <- 1./(se.x * se.x)
        w2.y     <- 1./(se.y * se.y)
        invsumw2 <- 1./(w2.x + w2.y)
        mbeta    <- (b.x*w2.x + b.y*w2.y)*invsumw2
        mse      <- sqrt(invsumw2)
        mp       <- pchisq(mbeta*mbeta/(mse*mse), 1, lower.tail=FALSE)

        if (name.x != "POOLED") {
            npops <- 1*(!is.na(data$beta.x)) + 1*(!is.na(data$beta.y))
        } else {
            npops <- data$npops + 1*(!is.na(data$beta.y))
            if (any(is.na(npops))) {
                npops[is.na(npops)] <-
                    1 * (!is.na(data$beta.x[is.na(npops)])) +
                        1 * (!is.na(data$beta.y[is.na(npops)]))
            }
        }
        na.x <- which(is.na(b.x) & !is.na(b.y))
        na.y <- which(is.na(b.y) & !is.na(b.x))
#	na.xy <- which(is.na(b.y) & is.na(b.x))
#	if (length(na.xy)>0) stop("Number of subjects missing in both data sets...")
        out           <- data.frame(name=data$name, stringsAsFactors=FALSE)
        out$strand    <- data$strand.x
        out$allele1   <- data$allele1.x
        out$allele2   <- data$allele2.x
        out$effallele <- data$effallele.x

        if (chr.x & chr.y) {
            out$chromosome <- data$chromosome.x
            out$chromosome[na.x] <- data$chromosome.y[na.x]
        } else if (chr.x) {
            out$chromosome <- data$chromosome.x
        } else if (chr.y) {
            out$chromosome <- data$chromosome.y
        }
#	print(c(pos.x, pos.y))

        if (pos.x & pos.y) {
            out$position <- data$position.x
            out$position[na.x] <- data$position.y[na.x]
        } else if (pos.x | pos.y) {
            out$position <- data$position
        }

        out$strand[na.x] <- data$strand.y[na.x]
        out$allele1[na.x] <- data$allele1.y[na.x]
        out$allele2[na.x] <- data$allele2.y[na.x]
        out$effallele[na.x] <- data$effallele.y[na.x]
        out$n <- data$n.x + data$n.y
        out$n[na.x] <- data$n.y[na.x]
        out$n[na.y] <- data$n.x[na.y]
        out$npops <- npops
        out$beta <- mbeta
        out$sebeta <- mse
        out$beta[na.x] <- b.y[na.x]
        out$sebeta[na.x] <- se.y[na.x]
        out$beta[na.y] <- b.x[na.y]
        out$sebeta[na.y] <- se.x[na.y]

        na.x <- which(is.na(data$effallelefreq.x))
        na.y <- which(is.na(data$effallelefreq.y))
        nNAb <- which(!is.na(data$effallelefreq.x) &
                      !is.na(data$effallelefreq.y))

        out$effallelefreq <- rep(NA, dim(data)[1])
        out$effallelefreq[nNAb] <- (data$effallelefreq.x[nNAb] *
                                    data$n.x[nNAb] +
                                    data$effallelefreq.y[nNAb] *
                                    data$n.y[nNAb]) / (data$n.x[nNAb]
                                                       +
                                                       data$n.y[nNAb])
        out$effallelefreq[na.x] <- data$effallelefreq.y[na.x]
        out$effallelefreq[na.y] <- data$effallelefreq.x[na.y]
#	print(data$effallelefreq.x)
#	print(data$effallelefreq.y)
#	print(out$effallelefreq)

        na.x <- which(is.na(data$call.x))
        na.y <- which(is.na(data$call.y))
        nNAb <- which(!is.na(data$call.x) & !is.na(data$call.y))
        out$call <- rep(NA, dim(data)[1])
#	print(na.x)
#	print(nNAb)
        out$call[nNAb] <- (data$call.x[nNAb] * data$n.x[nNAb] +
                           data$call.y[nNAb] * data$n.y[nNAb]) /
                               (data$n.x[nNAb] + data$n.y[nNAb])
        out$call[na.x] <- data$call.y[na.x]
        out$call[na.y] <- data$call.x[na.y]

#	print(data$pexhwe.y)
        if (all(data$pexhwe.x <= 1.0, na.rm=TRUE)) {
            data$pexhwe.x <- qchisq(1. - data$pexhwe.x, 1)
        }
        if (all(data$pexhwe.y <= 1.0, na.rm=TRUE)) {
            data$pexhwe.y <- qchisq(1. - data$pexhwe.y, 1)
        }
#	print(data$pexhwe.y)
        na.x <- which(is.na(data$pexhwe.x))
        na.y <- which(is.na(data$pexhwe.y))
        nNAb <- which(!is.na(data$pexhwe.x) & !is.na(data$pexhwe.y))
        out$pexhwe <-  rep(NA, dim(data)[1])
        out$pexhwe[nNAb] <- data$pexhwe.x[nNAb]+data$pexhwe.y[nNAb]
        out$pexhwe[na.x] <- data$pexhwe.y[na.x]
        out$pexhwe[na.y] <- data$pexhwe.x[na.y]

        s <- match(substr(names(data), 1, 5), "obeta")
        if (any(!is.na(s))) s <- which(!is.na(s))
        if (any(!is.na(s))) out <- cbind(out, data[, names(data)[s]])
        if (name.x != "POOLED") out[, paste("obeta", name.x,sep="")] <- b.x

        out[, paste("obeta", name.y,sep="")] <- b.y
        s <- match(substr(names(data), 1, 3), "ose")

        if (any(!is.na(s))) s <- which(!is.na(s))
        if (any(!is.na(s))) out <- cbind(out, data[, names(data)[s]])
        if (name.x != "POOLED") out[, paste("ose", name.x,sep="")] <- se.x

        out[, paste("ose", name.y,sep="")] <- se.y
        cat("Lambda POOLED data =", median((out$beta / out$sebeta)^2,
                                           na.rm=TRUE) / chi2med, "\n")
        out$chi2 <- out$beta * out$beta / (out$sebeta * out$sebeta)
        out$p <- pchisq(out$chi2, 1, lower.tail=FALSE)
        out
}


swap.y <- function(df, swcode) {
        df$beta.y[swcode] <- (-1.) * df$beta.y[swcode]
        df$effallelefreq.y[swcode] <- (1. - df$effallelefreq.y[swcode])
        a1.save <- df$allele1.y[swcode]
        a2.save <- df$allele2.y[swcode]
        df$effallele.y[swcode] <- a1.save
        df$allele1.y[swcode]   <- a2.save
        df$allele2.y[swcode]   <- a1.save
        df
}


swap.x <- function(df, swcode) {
        df$beta.x[swcode] <- (-1.) * df$beta.x[swcode]
        df$effallelefreq.x[swcode] <- (1. - df$effallelefreq.x[swcode])
        a1.save <- df$allele1.x[swcode]
        a2.save <- df$allele2.x[swcode]
        df$effallele.x[swcode] <- a1.save
        df$allele1.x[swcode]   <- a2.save
        df$allele2.x[swcode]   <- a1.save
        df
}


recode.y <- function(df, swcode) {
        a2t <- swcode[df$effallele.y[swcode]=="A"]
        t2a <- swcode[df$effallele.y[swcode]=="T"]
        g2c <- swcode[df$effallele.y[swcode]=="G"]
        c2g <- swcode[df$effallele.y[swcode]=="C"]
        df$effallele.y[a2t] <-"T"
        df$effallele.y[t2a] <-"A"
        df$effallele.y[g2c] <-"C"
        df$effallele.y[c2g] <-"G"
        a2t <- swcode[df$allele1.y[swcode]=="A"]
        t2a <- swcode[df$allele1.y[swcode]=="T"]
        g2c <- swcode[df$allele1.y[swcode]=="G"]
        c2g <- swcode[df$allele1.y[swcode]=="C"]
        df$allele1.y[a2t] <-"T"
        df$allele1.y[t2a] <-"A"
        df$allele1.y[g2c] <-"C"
        df$allele1.y[c2g] <-"G"
        a2t <- swcode[df$allele2.y[swcode]=="A"]
        t2a <- swcode[df$allele2.y[swcode]=="T"]
        g2c <- swcode[df$allele2.y[swcode]=="G"]
        c2g <- swcode[df$allele2.y[swcode]=="C"]
        df$allele2.y[a2t] <-"T"
        df$allele2.y[t2a] <-"A"
        df$allele2.y[g2c] <-"C"
        df$allele2.y[c2g] <-"G"
        f2r <- swcode[df$strand.y[swcode]=="+"]
        r2f <- swcode[df$strand.y[swcode]=="-"]
        t2b <- swcode[df$strand.y[swcode]=="TOP"]
        b2t <- swcode[df$strand.y[swcode]=="BOT"]
        df$strand.y[f2r] <- "-"
        df$strand.y[r2f] <- "+"
        df$strand.y[t2b] <- "BOT"
        df$strand.y[b2t] <- "TOP"

        if (any(is.na(df$strand.y[swcode]))) {
            warning("NAs in strand!")
        }

        if (any(df$strand.y[swcode]!="+" &
                df$strand.y[swcode]!="-" &
                df$strand.y[swcode]!="TOP" &
                df$strand.y[swcode]!="BOT")) {
            warning("Unrecognised coding in strand!")
        }
        df
}
