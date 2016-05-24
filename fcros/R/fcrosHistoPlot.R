fcrosHistoPlot <- function(af, nbins = 50) {
    hist(af$ri, nclass = nbins, xlab = "", main = "FCROS statistics", xlim = c(0, 1));
}
