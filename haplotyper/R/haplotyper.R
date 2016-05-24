#' haplotyper function identifies haplotypes within QTL.
#'
#' This function groups together all individuals of a population
#' with the same haplotype.
#'
#' @param x  a data.frame that should be loaded with read.table function.
#' Each row represents the individuals while each column represents the markers.
#' The first column contains the names of the genotypes.
#' @param Print option for print the haplotyper result.
#' The default is FALSE
#'
#' @return a matrix with the haplotypes
#' @author Sebastian Simondi, Victoria Bonnecarrere, Lucia Gutierrez,
#' Gaston Quero
#' @details Each group contains individual with the same allele in each SNP,
#'whether or not missing data.
#' @seealso read.table function
#'
#' @export
#'
#' @import graphics
#'
#' @import utils
#'
#' @examples
#' \dontrun{
#' data(rice_qtl)
#' haplotyper(rice_qtl)
#'}
haplotyper <- function(x, Print = FALSE) {
    H_data_o <- x
    H_data <- x
    dir.create("haplotyper_reports", showWarnings = F)
    hplq.result <- NULL
    hplq.final <- NULL
    id <- H_data[, 1]
    Q <- H_data[, -1]
    hplq <- NULL
    hp <- NULL
    hp.1 <- NULL
    cq <- NULL
    w <- NULL
    cr <- NULL
    c.Q <- rowSums(Q == 0)
    y <- which(c.Q == 0)
    b <- Q[y, ]
    Y <- b[!duplicated(b), ]
    for (i in 1:nrow(Q)) {
        for (j in 1:nrow(Y)) {
            cq <- rowSums(Q[i, ] == 0)
            w <- Q[i, ] - Y[j, ]
            cr <- rowSums(w == 0)
            zeros <- cq + cr
            if (zeros == ncol(Q)) {
                hp <- cbind(hp, i)
                hp.1 <- cbind(hp.1, j)
                hpl <- cbind(Q[i, ], j)
                hpl1 <- cbind(id[i], hpl)
                hplq <- rbind(hplq, hpl1)
            }
        }
    }
    u <- H_data[-c(hp), ]
    fa <- NULL
    for (i in 1:nrow(Y)) {
        fa <- (cbind(fa, rowSums(hp.1 == i)))
    }
    fa <- cbind(fa, nrow(u))
    fr <- (fa / nrow(Q)) * 100
    fh <- round(rbind(fa, fr), 1)
    colnames(fh) <- c(paste("haplo", 1:nrow(Y)), "undetermined")
    freq <- c("freq.abs", "freq.rel")
    fh <- data.frame(freq, fh)
    hap.id <- paste("haplo", 1:nrow(Y))
    Y <- cbind(hap.id, Y)
    colnames(hplq)[colnames(hplq) == "id[i]"] <- "id.geno"
    colnames(hplq)[colnames(hplq) == "j"] <- "haplo.qtl"
    undetermined <- cbind(u, "undetermined")
    colnames(undetermined) <- colnames(hplq)
    hplq.final <- rbind(hplq, undetermined)

    d <- duplicated (hplq.final[, 1]) | duplicated(hplq.final[, 1],
                    fromLast = TRUE)
    d <- hplq.final[d, ]

    hplq.result <- list(h.result = hplq.final, haplotypes = Y,
                        duplicates = d, freq = fh, und = u)

    write.table(hplq.result$freq,
                file = "./haplotyper_reports/haplotypes_frequencies.txt",
                sep = ",", eol = "\n", na = "-", dec = ".",
                col.names = T, row.names = F)

    write.table(hplq.result$haplotypes,
                file = "./haplotyper_reports/haplotypes.txt",
                sep = ",", eol = "\n", na = "-", dec = ".",
                col.names = T, row.names = F)

    write.table(hplq.result$u,
                file = "./haplotyper_reports/undeterminates.txt", sep = ",",
                eol = "\n", na = "-", dec = ".",
                col.names = T, row.names = F)

    write.table(hplq.result$duplicates,
                file = "./haplotyper_reports/duplicates.txt",
                sep = ",", eol = "\n", na = "-", dec = ".",
                col.names = T, row.names = F)

    write.table(hplq.result$h.result,
                file = "./haplotyper_reports/haplotyper_result.txt",
                sep = ",", eol = "\n", na = "-", dec = ".",
                col.names = T, row.names = F)

    haplo.image <- hplq.final[order(hplq.final$haplo.qtl, decreasing = TRUE), ]
    haplo.image$id.geno <- NULL
    haplo.image$haplo.qtl <- NULL
    haplo.image[haplo.image == 0] <- NA

    layout(mat = matrix(c(1, 2), 1, 2, byrow = TRUE))
    input <- H_data_o[, -1]
    input[input == 0] <- NA

    image(1:ncol(input), 1:nrow(input), t(input),
          ylab = "Individuals", xlab = "Markers",
          col = c("red2", "orange", "darkgreen", "navyblue"), main = "input")
    box()

    image(1:ncol(haplo.image), 1:nrow(haplo.image), t(haplo.image),
          ylab = "Individuals",xlab = "Markers",
          col = c("red2", "orange", "darkgreen", "navyblue"), main = "output")
    box()
    op <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
              mar = c(0, 0, 0, 0), new = TRUE)

    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

    legend("bottom", legend = c("A", "C", "G", "T"), cex = 0.9,
           bty = "n", text.col = "black",
           pch = 15, pt.cex = 2,
           col = c("red2", "orange", "darkgreen", "navyblue"),
           pt.lwd = 2,horiz = TRUE)
    par(op)
    if (Print == TRUE) {
        print(hplq.result)
    }

}
