`plot.WGassociation` <-
function (x, alpha = 0.05, plot.all.SNPs = FALSE, print.label.SNPs = TRUE, 
    cutPval = c(0, 1e-10, 1), whole, ylim.sup = ifelse(is.null(attr(x,"fast")),1e-40,
    1e-30), col.legend = c("red","gray60"), sort.chromosome=TRUE, centromere, ...) 
{
    if (!inherits(x, "WGassociation")) 
        stop("x must be an object of class 'WGassociation'")

    if (missing(whole))
     {
       if (ncol(attr(x, "pvalues"))==2 & length(unique(attr(x,"gen.info")[,2]))>10)
         whole <- TRUE
       else 
         whole <- FALSE 
     }
    if (!whole) {
        x <- attr(x, "pvalues")
        ylims<-range(-log(x[,-1],10),na.rm=TRUE)
        control.Mono <- grep("Monomorphic", as.character(x[, 
            1]))
        control.Geno <- grep("Genot", as.character(x[, 1]))

        if (plot.all.SNPs | length(c(control.Mono, control.Geno))==0) 
            ans <- x
        else ans <- x[-c(control.Mono, control.Geno), ]
        old.mar <- par("mar")
        old.mfrow <- par("mfrow")
        on.exit(par(mar = old.mar, mfrow = old.mfrow))
#VM
        n.models <- ncol(ans) - 1
#        m <- matrix(c(1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 8, 1, 9, 1, 10, 1, 11, 1, 12), nrow = 11, ncol = 2, byrow = TRUE)        

        m <-cbind(rep(1,n.models*2+1),2:(n.models*2+2))
                        
        layout(m, heights = c(rep(c(0.3, 1), n.models), 0.7), widths = c(0.05, 1))
# para que salgan bien los labels de los SNPs si no se dibujan todos los modelos       

        par(mar = c(0, 0, 0, 0))
# region vertical
#        if (n.models>2)
#         control.labelY <- 6.5 - ceiling(n.models/2)
#        else if (n.models==2)
#         control.labelY <- 4.5
#        else
#          control.labelY <- 5

        control.labelY <- 3.

        plot(rep(1, 5), 1:5, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, control.labelY, expression(paste(-log[10]," (p value)")), font = 2, srt = 90, cex = 2, adj=0.5)
#        models <- c("codominant", "dominant", "recessive", "overdominant", "log-additive")
        models <- names(ans)[-1]
        ok <- 2
        for (i in 1:n.models) {
            par(mar = c(0, 0, 0, 0))
         #   if (any(models[i] %in% names(ans))) {
            # cabecera
                plot(1:3, rep(1, 3), type = "n", axes = FALSE, 
                  xlab = "", ylab = "")
                if (i == 1) 
                  legend(2.5, 1.5, c("Nominal p value", "Bonferroni correction"), 
                    lty = c(2, 2), col = c("pink1", "red"), cex = 1, 
                    bty = "n", y.intersp = 0.8)
                text(2, 1, models[i], col = "blue", cex = 1.5, 
                  font = 2)
                pval <- ans[, ok]
                xx <- -log(pval,10)

                par(mar = c(0.3, 0, 0, 0))
             # p-values   
                plot(c(1:length(xx)), xx, type = "b", axes = FALSE, 
                  ylab = "", xlab = "", ylim=ylims)
                axis(1, at = c(1:length(xx)), labels = rep("", 
                  length(xx)), pos = 0)
                control.y <- ceiling(seq(0,ceiling(ylims[2]),length=6))


                axis(2, at = control.y, labels = -control.y, 
                  pos = 1,xaxt = "n")

           
                cut <- -log(alpha,10)
                segments(1, cut, length(xx), cut, col = "pink2", 
                  lty = 2)
                pvalues <- xx[!is.na(xx)]
                cut.p <- alpha/length(pvalues)
                cut.trans <- -log(cut.p,10)
                if (cut.trans > max(xx, na.rm = TRUE)) {
                  cat("Warning: No SNP is statistically significant after \n         Bonferroni Correction under",
                    models[i], "model \n")
                }
#                else {
                  segments(1, cut.trans, length(xx), cut.trans, 
                    col = "red", lty = 2)
#                }
                ok <- ok + 1
        #    }
        }
        if (ok < n.models) 
            plot(1:3, rep(1, 3), type = "n", axes = FALSE, xlab = "", 
                ylab = "")
        plot(c(1:nrow(ans)), rep(0, nrow(ans)), type = "n", axes = FALSE)
        if (print.label.SNPs) {
            text(c(1:nrow(ans)), rep(1, nrow(ans)), dimnames(ans)[[1]], 
                srt = 90, adj = 1, xpd = TRUE, ...)
        }
        else {
            text(nrow(ans)/2, 0, "SNPs", cex = 2)
        }
    }
    else {
        gen.info <- attr(x, "gen.info")
        chrs <- gen.info[, 2]
        chr.l <- unique(chrs)

        if (sort.chromosome)
          chr <- chr.l[orderChromosome(chr.l)]
        else
         chr<-unique(chrs)

        pval <- attr(x, "pvalues")
        limits <- c(min(gen.info[, 3]), max(gen.info[, 3]))
        if (missing(centromere)) {
            centro <- c(1.23e+08, 93500000, 91500000, 5.1e+07, 
                47700000, 60500000, 58900000, 4.5e+07, 4.9e+07, 
                4e+07, 5.3e+07, 35400000, 1.5e+07, 15600000, 
                1.7e+07, 3.9e+07, 2.4e+07, 1.6e+07, 28500000, 
                27800000, 12200000, 11900000, 58500000, 1e+07)
        }
        else {
            centro <- centromere
        }
        n.chr <- length(chr)
        old.mfrow <- par("mfrow")
        old.mar <- par("mar")
        old.xpd <- par("xpd")
        on.exit(par(mfrow = old.mfrow, mar = old.mar))
        par(mfrow = c(26, 1))
        par(mar = c(0.5, 5, 0, 3))
        par(xpd = TRUE)
        plot(c(1:3), rep(1, 3), axes = FALSE, xlab = "", ylab = "", 
            type = "n")
        text(2, 1, paste("Genetic model:", names(pval)[2]), adj = 0.5, 
            font = 2, cex = 2)
        plot(c(1:3), rep(1, 2, 3), axes = FALSE, xlab = "", ylab = "", 
            type = "n")
        text(1, 1.2, "p value", adj = 0)
        legend(1.2, 1.2, levels(cut(runif(100), cutPval)), col = col.legend, 
            pt.bg = col.legend, pch = rep(22, 4), horiz = TRUE, 
            cex = 1, bty = "n", pt.cex = 1.6, yjust = 0.5)
        max.y <- min(pval[pval[,2]>0, 2], na.rm = TRUE)
        if (ylim.sup > max.y) 
            ylim.sup <- max.y
        for (i in 1:n.chr) {
            snp <- gen.info[gen.info[, 2] == chr[i], 1]
            pos <- gen.info[gen.info[, 2] == chr[i], 3]
            dat <- pval[dimnames(pval)[[1]] %in% snp, -1]
            col <- as.character(cut(dat, cutPval, labels = col.legend))
            col[is.na(col)] <- "white"
            dat[is.na(dat)] <- 1
            plot(pos, -log(dat,10), axes = FALSE, xlab = "", ylab = "", 
                type = "h", col = col, xlim = limits, ylim = c(0, 
                  -log(ylim.sup,10)))
            segments(limits[1], 0, max(pos), 0)
            segments(centro[i], -1, centro[i], -log(max.y,10), lwd = 3, 
                col = "darkblue", xpd = TRUE)
            text(par("usr")[1], 0, chr[i], cex = 1, adj = 0.5)
        }
        plot(limits, c(1, 1), axes = FALSE, xlab = "", ylab = "", 
            type = "n", col = col, xlim = limits)
        text(limits[1], 1, limits[1], cex = 1, adj = 0)
        text(limits[2], 1, limits[2], cex = 1, adj = 1)
        text((limits[2] - limits[1])/2, 1, "Genomic Position", 
            cex = 1, adj = 0, font = 2)
    }
}

