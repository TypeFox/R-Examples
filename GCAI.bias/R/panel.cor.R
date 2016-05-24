panel.cor <-
function(x, y, digits=2, prefix="", cex.cor, cor.cut=1)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x[which(x>=cor.cut|y>=cor.cut)], y[which(x>=cor.cut|y>=cor.cut)],use="pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * r)
}
