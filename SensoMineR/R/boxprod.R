boxprod<-function(donnee,col.p,firstvar,lastvar=ncol(donnee),numr = 2,numc = 2) {

nbdesc <- lastvar-firstvar+1

mult <- nbdesc %/% (numr*numc)

if (nbdesc==(nbdesc %/% (numr*numc))*(numr*numc)) mult=mult-1
for (m in 0:mult) {
    par(mfrow = c(numr,numc))

    for (nbd in 1:(numr*numc)) {
    nb <- (m*(numr*numc)+nbd)
          if (nb <= nbdesc)       {
boxplot(donnee[,nb+firstvar-1]~donnee[,col.p],col="orchid3",main = names(donnee[nb+firstvar-1]),boxfill= "light gray",outpch = 21:25, outlty = 2,bg = "pink",lwd = 0.5, medcol = "dark blue", medcex = 1, medpch=15)
                              } 
                   }
if (m < mult) dev.new()
          }  #for (m in 0:mult) {

                                  }
