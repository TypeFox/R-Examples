hist.kselect <- function(x, xax = 1, mar=c(0.1,0.1,0.1,0.1),
                         ncell=TRUE, csub=2,
                         possub=c("bottomleft", "topleft",
                         "bottomright", "topright"),
                         ncla=15, ...)
{
    ## Verifications
    possub<-match.arg(possub)
    if (!inherits(x, "kselect")) stop("should be a 'kselect' object")


    ## 1. Creation de la liste
    Xi<-x$initab
    Xrecalc<-t(as.matrix(apply(Xi, 1,
                               function(y) y*x$lw/sum(x$lw))))%*%as.matrix(x$l1)
    li.Xi<-split(as.data.frame(Xrecalc), x$initfac)
    li.wei<-split(x$initwei, x$initfac)
    rx<-range(Xrecalc[,xax])
    br<-seq(rx[1]-(rx[2]-rx[1])/100, rx[2]+(rx[2]-rx[1])/100, length=ncla)

    ## Graphical settings
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    ngraph<-length(li.Xi)
    par(mfrow = n2mfrow(ngraph), mar=mar)

    lapply(1:ngraph, function(i) {
        Xtmp<-li.Xi[[i]][,xax]
        poids<-li.wei[[i]]
        poids[poids>0]<-1
        histniche(data.frame(Xtmp), poids, axes=FALSE, ncla=ncla,
                  main="",...)
        scatterutil.sub(names(li.Xi)[i],
                        csub=csub, possub=possub)
        box()
    })


}
