panel.cor <- function(x, y, digits = 2, cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y,use="na.or.complete")
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}

