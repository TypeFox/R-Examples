abuocc <- function (taxa, minabu = 0, panel='all')
{
    if (!is.data.frame(taxa))
        taxa <- data.frame(taxa)
    spc.plt <- apply(taxa > minabu, 1, sum)
    plt.spc <- apply(taxa > minabu, 2, sum)
    if (minabu == 0) {
        mean.abu <- apply(taxa, 2, sum)/plt.spc
    } else {
        mean.abu <- rep(0, ncol(taxa))
        for (i in 1:ncol(taxa)) {
            mask <- taxa[, i] > minabu
            mean.abu[i] <- sum(taxa[mask, i])/max(1, plt.spc[i])
        }
    }
    mean.abu[is.na(mean.abu)] <- 0
    if (panel=='all' || panel==1) {
        plot(rev(sort(plt.spc[plt.spc > minabu])), log = "y", xlab = "Species Rank",
            ylab = "Number of Plots", main = "Species Occurrence")
        if (panel == 'all') readline("Press return for next plot ")
    }
    if (panel=='all' || panel==2) {
        plot(rev(sort(spc.plt)), xlab = "Plot Rank", ylab = "Number of Species",
            main = "Species/Plot")
        if (panel=='all') readline("Press return for next plot ")
    }
    if (panel=='all' || panel==3) {
        plot(plt.spc[mean.abu > minabu], mean.abu[mean.abu > minabu],
            log = "y", xlab = "Number of Plots", ylab = "Mean Abundance",
            main = "Abundance vs Occurrence")
        yorn <- readline("Do you want to identify individual species? Y/N : ")
        if (yorn == "Y" || yorn == "y")
            identify(plt.spc[mean.abu > minabu], mean.abu[mean.abu >
                minabu], names(taxa)[mean.abu > minabu])
        if (panel=='all') readline("Press return for next plot ")
    }
    if (panel=='all' || panel==4) {
        plot(spc.plt, apply(taxa, 1, sum), xlab = "Number of Species/Plot",
            ylab = "Total Abundance")
        yorn <- readline("Do you want to identify individual plots? Y/N : ")
        if (yorn == "Y" || yorn == "y")
            identify(spc.plt, apply(taxa, 1, sum), labels = row.names(taxa))
    }
    out <- list(spc.plt = spc.plt, plt.spc = plt.spc, mean = mean.abu)
    attr(out,'call') <- match.call()
    attr(out,'taxa') <- deparse(substitute(taxa))
    invisible(out)
}
