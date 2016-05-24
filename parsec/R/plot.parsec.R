plot.parsec <-
function(x, which = c("Hasse", "threshold", "identification", "rank", "gap"), ask = dev.interactive(),  shape = c("square", "circle", "equispaced"), ...) {
    
    if(!all(is.numeric("labels"), is.character("labels")))
        labels <- rownames(x$incidence)
    
    avgrank <- apply(
        x$rank_dist,
        2,
        function(v) sum((1:x$number_of_profiles)*v)
    )
    ord <- order(x$idn_f, -avgrank)
        
    # diagramma di Hasse
    if (any(c("Hasse", "hasse") %in% which)) {
    plot(
        x$incidence,
        shape = shape,
        pch = 21,
        cex = max(nchar(rownames(x$incidence)))+2,
        bg = "white",
        lwd = 1 + 2*x$threshold
    )
    
    if(ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    }
    
    # frequenze relative della threshold
    # per valutarne l'influenza sui risultati
    if ("threshold" %in% which) {
    plot(
        x$thr_dist[x$threshold],
        type = "h",
        axes = FALSE, frame.plot = TRUE,
        xlab = "Threshold",
        ylab = "Relative frequencies",
        xlim = c(0, sum(x$threshold) + 1),
        ylim = c(0, 1),
        lwd = 3
    )
    axis(1, at = 1:sum(x$threshold), labels = labels[x$threshold])
    axis(2, at = 0:4/4)
        
    if(ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    }
    
    # identification function
    if ("identification" %in% which) {
    plot(
        x$idn_f[ord],
#         type = "h",
        type = "l",
        axes = FALSE,
        frame.plot = TRUE,
        ylim = 0:1,
        lwd = 3,
        xlab = "Profiles",
        ylab = "Identification function",
#         panel.first = {
#             lines(1:x$number_of_profiles, x$idn_f[ord], lty = 2, col = "gray")
#         }
    )
    axis(1, at = 1:x$number_of_profiles, labels = labels[ord])
    axis(2, at = 0:4/4)
    
    if(ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    }
    
    # frequenze relative dei ranghi
    if ("rank" %in% which) {
    barplot(
        t(x$rank_dist[ord,]),
        col=gray(1-0:x$number_of_profiles/x$number_of_profiles),
        xlab = "Profiles",
        ylab = "Rank distribution"
        # qua forse si dovrebbe stare attenti alle lables
    )
    
    if(ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    }
    
    try(
    {
    # gap relativi
    if ("gap" %in% which) {
    plot(-x$point_relative_poverty_gap[ord],
         type="h",
         ylim=c(-1, 1),
         axes = FALSE,
         frame.plot = TRUE,
         col = "black",
         lwd = 3,
         ylab ="Relative gap",
         xlab = "Profiles",
         panel.first={
             abline(v=1:x$number_of_profiles, lty = 2, col=ifelse(x$threshold[ord], "black", "gray"), lwd = 1)
             abline(h=c(0, -x$poverty_gap, x$wealth_gap), lty = c(1, 2, 2), col = c("black", "gray", "gray"))
         })
    points(
        x$point_relative_wealth_gap[ord],
        type="h",
        col = "black",
        lwd = 3
    )
    axis(1, 1:x$number_of_profiles, labels = labels[ord])
    axis(2, at=-2:2/2, labels = abs(-2:2/2))
    }
    }, silent = TRUE)
}
