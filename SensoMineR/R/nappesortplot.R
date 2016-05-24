#représentation graphique des nappes catégorisées : fonction nappecatplot
nappesortplot=function (donnee,first="nappe", numr = 2, numc = 2, lim = c(60, 40)) {
    nbjuge <- ncol(donnee)/3
#Si pour chaque sujet d'abord la catégorisation  puis la nappe
if (first=="catego"){
  donnee2=donnee
  for (i in 1:nbjuge){
  donnee2[,(3*(i-1)+1):(3*(i-1)+2)]=donnee[,(3*(i-1)+2):(3*i)]
  donnee2[,(3*i)]=donnee[,(3*(i-1)+1)]
  colnames(donnee2)[(3*(i-1)+1):(3*(i-1)+2)]=colnames(donnee)[(3*(i-1)+2):(3*i)]
  colnames(donnee2)[(3*i)]=colnames(donnee)[(3*(i-1)+1)]
  }
  donnee=donnee2
}
 palette = palette(c("black", "red", "green3", "blue", "cyan", "magenta",
            "darkgray", "darkgoldenrod", "darkgreen", "violet",
            "turquoise", "orange", "lightpink", "lavender", "yellow",
            "lightgreen", "lightgrey", "lightblue", "darkkhaki",
            "darkmagenta", "darkolivegreen", "lightcyan", "darkorange",
            "darkorchid", "darkred", "darksalmon", "darkseagreen",
            "darkslateblue", "darkslategray", "darkslategrey",
            "darkturquoise", "darkviolet", "lightgray", "lightsalmon",
            "lightyellow", "maroon"))

    mult <- nbjuge%/%(numr * numc)
    if (nbjuge == (nbjuge%/%(numr * numc)) * (numr * numc))
        mult = mult - 1
    for (m in 0:mult) {
        par(mfrow = c(numr, numc))
        for (nbd in 1:(numr * numc)) {
            nb <- (m * (numr * numc) + nbd)
            if (nb <= nbjuge) {
                plot(donnee[, (3 * (nb - 1) + 1):(3 * (nb - 1) + 2)],
                  xlim = c(0, lim[1]), ylim = c(0, lim[2]), xlab = "",
                  ylab = "", main = paste(colnames(donnee)[3 *
                    nb]), type = "n", asp = 1)
for (i in 1:nlevels(donnee[,3*nb])){
numl=which(donnee[,3*nb]==levels(donnee[,3*nb])[i])
                points(donnee[numl, (3 * (nb - 1) + 1):(3 * (nb - 1) + 2)],
                  col = i, cex = 0.8, pch = 20)
                text(donnee[numl, (3 * (nb - 1) + 1):(3 * (nb - 1) + 2)], label = rownames(donnee)[numl],
                  col = i, cex = 0.8, pos = 3, offset = 0.2)

h=mean(donnee[numl,(3*(nb-1)+1)])
k=mean(donnee[numl,(3*(nb-1)+2)])
a=max(abs(max(donnee[numl,(3*(nb-1)+1)])-h),abs(h-min(donnee[numl,(3*(nb-1)+1)])))+2.5
b=max(abs(max(donnee[numl,(3*(nb-1)+2)])-k),abs(k-min(donnee[numl,(3*(nb-1)+2)])))+2.5

x=seq(from=(-a+h),to=(a+h),by=0.01)
#enlève le premier et dernier élément sinon quelques problèmes dans la racine carrée pour y
x=x[-c(1,length(x))]
y=b*sqrt(1-((x-h)/a)^2)+k
x2=rev(x)
y2=-b*sqrt(1-((x2-h)/a)^2)+k



#lines(x,y,col=i)
#lines(x2,y2,col=i)
}

}
        }
        if (m < mult)
            dev.new()
    }
    par(las = 0)
}