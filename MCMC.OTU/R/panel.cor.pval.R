panel.cor.pval <- function(x, y, digits = 2, cex.cor, p.cut=0.1)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor.test(x,y,use="na.or.complete")$p.value
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    if(r>p.cut) { txt=""}
    text(0.5, 0.5, txt)
}
