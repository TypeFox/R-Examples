calc.Rsq <- function (cross, peak, th = 0.001, round) 
{
    require(qtl)
    if (length(class(cross)) < 2 || class(cross)[2] != "cross") 
        stop("Input should have class \"cross\".")
    if (!all(attr(peak, "class", exact = TRUE) %in% c("peak", 
        "list"))) 
        stop("Input should have class \"peak\".")
    if (!is.numeric(th) & !is.vector(th) & (th > 1 | th < 0)) 
        stop("Input should be a numeric vector of length 1 with value between 0 and 1.")
    scanone <- get(attr(peak, "scanone", exact = TRUE))
    if (!all(levels(scanone$chr) == names(cross$geno))) 
        stop("Arguments peak and/or cross misspecified: attributes(peak)$scanone should describe a scanone object which has been performed on cross object")
    if (!all(scanone$pos == pseudo.map(cross))) 
        stop("Arguments peak and/or cross misspecified: attributes(peak)$scanone should describe a scanone object which has been performed on cross object")
    if (!missing(round) && !is.numeric(round)) 
        stop("Argument round misspecified: round should be an integer >= 0")
    res <- list(un = NA)
    all.rsq <- data.frame(qtl = NA, rsq = NA, pF = NA)
    for (i in seq(length(peak))) {
        trait <- names(peak[i])
        if (!names(peak[i]) %in% names(cross$pheno)) 
            if (names(peak[i]) == "lod") 
                stop("phenotype '", names(peak[i]), "' not found in cross object. You should rename the trait in peak object by the real trait name.")
            else stop(" phenotype '", names(peak[i]), "' not found in cross object.")
        col <- grep(paste("^", trait, "$", sep = ""), names(cross$phe))
        perf <- cross$phe[[col]]

        n <- 0
        data <- data.frame(perf = perf)

        resbytrait <- list(un = NA)
        for (y in 1:length(peak[[i]])) {
            if (!is.na(peak[[i]][y])) {
                chr <- names(peak[[i]][y])
                for (z in seq(length(as.vector(peak[[i]][[y]]$mname.peak)))) {
                  pos <- peak[[i]][[y]]$peak.cM[z]
                  marker <- find.flanking(cross, as.numeric(chr), 
                    as.numeric(pos))
                  geno <- pull.geno(cross)[, paste(marker$close)]
                  geno <- data.frame(geno)
                  n <- n + 1
                  attributes(geno)$names <- paste(trait, ".", 
                    y, ".", z, sep = "")
                  data <- cbind(data, geno)
                }
            }
        }

        if (ncol(data) < 2) {
            next
        }
        else {
            c <- 0
            for (y in seq(data$perf)) if (any(grep("NA", data[y, 
                ]))) 
                c <- c(c, y)

            if(!(length(c)==1&&c==0)) data <- data[-c, ]

		marker <- names(data[-1])
            write.lm <- function(marker.names) {
                col <- ""
                row <- ""
                for (i in seq(length(marker.names))) {
                  c <- combn(marker.names, length(marker.names[1:i]))
                  m <- ""
                  for (y in seq(ncol(c))) {
                    v <- paste(c[, y], collapse = ":")
                    m <- c(m, v)
                  }
                  row <- paste(m[-1], collapse = "+")
                  col <- c(col, row)
                }
                f <- paste(col[-1], collapse = "+")
                invisible(f)
            }

            modele <- write.lm(marker)

            f <- paste("data$perf~", modele, sep = "")

            a <- anova(lm(as.formula(f), data = data))
            Rsq <- function(a, th) {
                all.rsq <- ""
                qtl <- ""
                p <- ""
                CMR <- a$"Sum Sq"[nrow(a)]
                CMT <- sum(a$"Sum Sq")
                for (y in seq(nrow(a) - 1)) {
                  if (a$"Pr(>F)"[y] <= th) {
                    CME <- a$"Sum Sq"[y]
                    rsq <- CME/CMT
                    all.rsq <- c(all.rsq, rsq)
                    qtl <- c(qtl, row.names(a[y, ]))
                    p <- c(p, a$"Pr(>F)"[y])
                  }
                  else {
                    all.rsq <- c(all.rsq, NA)
                    qtl <- c(qtl, row.names(a[y, ]))
                    p <- c(p, NA)
                  }
                }
                invisible(data.frame(qtl = qtl[-1], rsq = all.rsq[-1], 
                  pF = (p[-1])))
            }
            rsq <- Rsq(a = a, th = th)
            if (!missing(round)) {
                rsq$rsq <- round(as.numeric(as.vector(rsq$rsq)), 
                  round)
                rsq$pF <- signif(as.numeric(as.vector(rsq$pF)), 
                  round)
            }
            all.rsq <- rbind(all.rsq, rsq)
        }
    }
    all.rsq <- all.rsq[-1, ]
    attributes(all.rsq)$class <- c("rsq", "data.frame")
    return(all.rsq)
}
